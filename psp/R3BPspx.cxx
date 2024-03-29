/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BPspx.h"
#include "R3BTGeoPar.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BGeoPspx.h"
#include "R3BMCStack.h"
#include "R3BPspxPoint.h"
#include "TClonesArray.h"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoCone.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPara.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <stdlib.h>

R3BPspx::R3BPspx()
    : R3BPspx("")
{
}

R3BPspx::R3BPspx(const TString& geoFile,
               const TGeoTranslation& trans,
               const TGeoRotation& rot,
               const Float_t z1,
               const Float_t z2,
               const Float_t z3)
    : R3BPspx(geoFile, { trans, rot }, z1, z2, z3)
{
}

R3BPspx::R3BPspx(const TString& geoFile,
               const TGeoCombiTrans& combi,
               const Float_t z1,
               const Float_t z2,
               const Float_t z3)
    : R3BDetector("R3BPspx", kPSP, geoFile, combi)
    , fZ1(z1)
    , fZ2(z2)
    , fZ3(z3)
    , fPspxCollection(new TClonesArray("R3BPspxPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BPspx::~R3BPspx()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fPspxCollection)
    {
        fPspxCollection->Delete();
        delete fPspxCollection;
    }
}

void R3BPspx::Initialize()
{
    FairDetector::Initialize();

    LOG(INFO) << "R3BPspx: initialisation";
    LOG(DEBUG) << "R3BPspx: Vol. (McId) " << gMC->VolId("PSP1Log");
}

void R3BPspx::SetSpecialPhysicsCuts()
{
    LOG(INFO) << "-I- R3BPspx: Adding customized Physics cut ... ";

    if (gGeoManager)
    {
        TGeoMedium* pSi = gGeoManager->GetMedium("silicon");
        if (pSi && 0 == 1)
        {
            // Setting processes for Si only
            gMC->Gstpar(pSi->GetId(), "LOSS", 3);
            gMC->Gstpar(pSi->GetId(), "STRA", 1.0);
            gMC->Gstpar(pSi->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pSi->GetId(), "COMP", 1.0);
            gMC->Gstpar(pSi->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pSi->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pSi->GetId(), "BREM", 1.0);
            gMC->Gstpar(pSi->GetId(), "HADR", 1.0);
            gMC->Gstpar(pSi->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "MULS", 1.0);
            gMC->Gstpar(pSi->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Si Only
            Double_t cutE = fCutE; // GeV-> 1 keV

            LOG(INFO) << "-I- R3BPspx: silicon Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Si
            gMC->Gstpar(pSi->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
    } //! gGeoManager
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BPspx::ProcessHits(FairVolume* vol)
{
    // 2 Simple Det PLane
    // get Info from DCH planes
    Int_t copyNo = -1;
    Int_t planeNr = -1;
    // Get the Geo info from MC Point
    gMC->CurrentVolID(copyNo);
    gMC->CurrentVolOffID(1, planeNr);

    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        // fTime   = gMC->TrackTime() * 1.0e09;
        // fLength = gMC->TrackLength();
        fTime_in = gMC->TrackTime() * 1.0e09;
        fLength_in = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Set additional parameters at exit of active volume. Create R3BPspxPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        if (fELoss == 0.)
            return kFALSE;

        fTime_out = gMC->TrackTime() * 1.0e09; // also in case particle is stopped in detector, or decays...
        fLength_out = gMC->TrackLength();
        fTime = (fTime_out + fTime_in) / 2.;
        fLength = (fLength_out + fLength_in) / 2.;

        if (gMC->IsTrackExiting())
        {
            const Double_t* oldpos;
            const Double_t* olddirection;
            Double_t newpos[3];
            Double_t newdirection[3];
            Double_t safety;

            gGeoManager->FindNode(fPosOut.X(), fPosOut.Y(), fPosOut.Z());
            oldpos = gGeoManager->GetCurrentPoint();
            olddirection = gGeoManager->GetCurrentDirection();

            for (Int_t i = 0; i < 3; i++)
            {
                newdirection[i] = -1 * olddirection[i];
            }

            gGeoManager->SetCurrentDirection(newdirection);
            // TGeoNode *bla = gGeoManager->FindNextBoundary(2);
            safety = gGeoManager->GetSafeDistance();

            gGeoManager->SetCurrentDirection(-newdirection[0], -newdirection[1], -newdirection[2]);

            for (Int_t i = 0; i < 3; i++)
            {
                newpos[i] = oldpos[i] - (3 * safety * olddirection[i]);
            }

            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
        }

        AddHit(fTrackID,
               fVolumeID,
               planeNr,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of PspxPoints for this track
        R3BStack* stack = (R3BStack*)gMC->GetStack();
        stack->AddPoint(kPSP);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BPspx::BeginEvent() {}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BPspx::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fPspxCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BPspx::Register() { FairRootManager::Instance()->Register("PSPPoint", GetName(), fPspxCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BPspx::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fPspxCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BPspx::Print(Option_t* option) const
{
    Int_t nHits = fPspxCollection->GetEntriesFast();
    LOG(INFO) << "R3BPspx: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BPspx::Reset()
{
    fPspxCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BPspx::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(INFO) << "R3BPspx: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BPspxPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = (R3BPspxPoint*)cl1->At(i);
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BPspxPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(INFO) << "R3BPspx: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BPspxPoint* R3BPspx::AddHit(Int_t trackID,
                            Int_t detID,
                            Int_t plane,
                            TVector3 posIn,
                            TVector3 posOut,
                            TVector3 momIn,
                            TVector3 momOut,
                            Double_t time,
                            Double_t length,
                            Double_t eLoss)
{
    TClonesArray& clref = *fPspxCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
    {
        LOG(INFO) << "R3BPspx: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    }
    return new (clref[size]) R3BPspxPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

// -----   Public method ConstructGeometry   ----------------------------------
void R3BPspx::ConstructGeometry()
{
    TString fileName = GetGeometryFileName();
    if (fileName.EndsWith(".root"))
    {
        LOG(INFO) << "Constructing PSP geometry from ROOT file " << fileName.Data();
        ConstructRootGeometry();

        // TODO: Check if this works as expected
        if (!fCombiTrans.IsIdentity())
        {
            TGeoNode* n = gGeoManager->GetTopNode()->GetDaughter(gGeoManager->GetTopNode()->GetNdaughters() - 1);
            TGeoCombiTrans* combtrans = (TGeoCombiTrans*)((TGeoNodeMatrix*)n)->GetMatrix();

            *combtrans = fCombiTrans;
        }

        TGeoNode* psp_node = gGeoManager->GetTopVolume()->GetNode("PSPWorld_0");

        TGeoNode* node = psp_node->GetVolume()->GetNode("PSP1LogWorld_1");
        TGeoCombiTrans* combtrans = (TGeoCombiTrans*)((TGeoNodeMatrix*)node)->GetMatrix();
        combtrans->SetDz(fZ1);

        node = psp_node->GetVolume()->GetNode("PSP2LogWorld_2");
        combtrans = (TGeoCombiTrans*)((TGeoNodeMatrix*)node)->GetMatrix();
        combtrans->SetDz(fZ2);

        node = psp_node->GetVolume()->GetNode("PSP3LogWorld_3");
        combtrans = (TGeoCombiTrans*)((TGeoNodeMatrix*)node)->GetMatrix();
        combtrans->SetDz(fZ3);
    }
    else
    {
        LOG(FATAL) << "PSP geometry file is not specified";
        exit(1);
    }
}

Bool_t R3BPspx::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("PSP1Log") || TString(name).Contains("PSP2Log") || TString(name).Contains("PSP3Log"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BPspx)
