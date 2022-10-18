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

#include "R3BPspxDigitizer.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

// includes for modeling
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoShapeAssembly.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>
#include <string>

#include "R3BMCTrack.h"
#include "R3BPspxPoint.h"

using std::cout;
using std::endl;

R3BPspxDigitizer::R3BPspxDigitizer()
    : FairTask("R3B Pspx Digitization scheme ")
{
}

R3BPspxDigitizer::~R3BPspxDigitizer() {}

void R3BPspxDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fPspxDigiPar = (R3BPspxDigiPar*)(rtdb->getContainer("R3BPspxDigiPar"));

    if (fPspxDigiPar)
    {
        cout << "-I- R3BPspxDigitizer::SetParContainers() " << endl;
        cout << "-I- Container R3BPspxDigiPar  loaded " << endl;
    }
}

InitStatus R3BPspxDigitizer::Init()
{

    //  cout<<"Init "<<endl;
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fPspxPoints = (TClonesArray*)ioman->GetObject("PSPPoint");
    fPspxMCTrack = (TClonesArray*)ioman->GetObject("MCTrack");

    // Register output array DchDigi
    fPspxDigi = new TClonesArray("R3BPspxDigi", 1000);
    ioman->Register("PspxDigi", "Digital response in Pspx", fPspxDigi, kTRUE);

    eventNoPspx = 0;

    // Initialise control histograms
    //      PspxXhis = new TH1F("PspxXhis","PspxXhis",600,-3.,3.);
    //      PspxXhis->GetXaxis()->SetTitle("Position");
    //      PspxXhis->GetYaxis()->SetTitle("Counts");

    return kSUCCESS;
}

void R3BPspxDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoPspx += 1;
    //     if(eventNoPspx/1000. == (int)eventNoPspx/1000.) cout<<"Event #: "<<eventNoPspx-1<<endl;

    Int_t nentriesPspx = fPspxPoints->GetEntries();
    Int_t TrackIdPspx = 0;

    Int_t psp3mul;
    Double_t psp3x, psp3y, psp3e;

    //******************** PSP3 ********************//
    psp3mul = 0;
    psp3x = 0.0;
    psp3y = 0.0;
    psp3e = 0.0;

    for (Int_t l = 0; l < nentriesPspx; l++)
    {
        //   cout<<"entries-Pspx "<<l<<endl;

        R3BPspxPoint* psp_obj = (R3BPspxPoint*)fPspxPoints->At(l);

        TrackIdPspx = psp_obj->GetTrackID();
        R3BMCTrack* aTrack = (R3BMCTrack*)fPspxMCTrack->At(TrackIdPspx);
        Int_t PID = aTrack->GetPdgCode();
        //     Int_t mother = aTrack->GetMotherId();

        Double_t fX_in = psp_obj->GetXIn();
        Double_t fY_in = psp_obj->GetYIn();
        Double_t fZ_in = psp_obj->GetZIn();
        Double_t fX_out = psp_obj->GetXOut();
        Double_t fY_out = psp_obj->GetYOut();
        Double_t fZ_out = psp_obj->GetZOut();
        Double_t PSPeloss = psp_obj->GetEnergyLoss() * 1000;

        Double_t fX = ((fX_in + fX_out) / 2);
        Double_t fY = ((fY_in + fY_out) / 2);
        Double_t fZ = ((fZ_in + fZ_out) / 2);

        //     if (PID==1000080150 && mother<0){
        if (PID > 1000501000)
        {
            if (fZ > 0.0)
            { // i.e. psp3

                //      psp3x=(((fX + 0.0) * 1.0) - ((fZ - 94.1) * (0.0)));
                psp3x = fX;
                psp3y = fY;
                psp3e += PSPeloss;

                //      PspxXhis->Fill(psp3x);
                psp3mul++;
            } // psp3
        }     // PID
    }

    // psp3x = gRandom->Gaus(psp3x, 0.0200);
    // psp3y = gRandom->Gaus(psp3y, 0.0200);

    AddHit(psp3mul, psp3x, psp3y, psp3e);
}
// -------------------------------------------------------------------------

void R3BPspxDigitizer::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fPspxDigi)
        fPspxDigi->Clear();
}

void R3BPspxDigitizer::Finish()
{
    // Write control histograms

    //    PspxXhis->Write();
}

R3BPspxDigi* R3BPspxDigitizer::AddHit(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e)
{
    TClonesArray& clref = *fPspxDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BPspxDigi(psp3mul, psp3x, psp3y, psp3e);
}

// R3BDchDigi* R3BDchDigitizer::AddHit(
// return new(clref[size]) R3BDchDigi();
//}

ClassImp(R3BPspxDigitizer)
