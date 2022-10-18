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

// ----------------------------------------------------------------
// -----                 R3BPspxMapped2CalPar              -----
// -----              Created  Jul 2022 by A. Jedele          -----
// -----             Adopted from script by S. Storck         -----
// -----  Purpose of script is to generate the 
// ----------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRuntimeDb.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxMapped2CalPar.h"
#include "R3BPspxMappedData.h"

#include "TClonesArray.h"
#include <iostream>
#include <limits>

R3BPspxMapped2CalPar::R3BPspxMapped2CalPar()
    : FairTask("PspxMapped2CalPar", 1)
    , fMapped(nullptr)
    , fCalItems(new TClonesArray("R3BPspxHitData"))
    , fCalPar(nullptr)
    , fTpat(-1)
    , fNumDet(0)
    , fNumStrips(0)
    , fUpdateRate(1000000)
    , fMinStats(100000)
{
}

R3BPspxMapped2CalPar::R3BPspxMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapped(nullptr)
    , fCalItems(new TClonesArray("R3BPspxHitData"))
    , fCalPar(nullptr)
    , fTpat(-1)
    , fNumDet(0)
    , fNumStrips(0)
    , fUpdateRate(1000000)
    , fMinStats(100000)
{
}

R3BPspxMapped2CalPar::~R3BPspxMapped2CalPar()
{
    delete fCalPar;
}

void R3BPspxMapped2CalPar::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan) {return kFATAL;}

    fMapped = (TClonesArray*)fMan->GetObject("PspxMapped");
    if (!fMapped) {return kFATAL;}

    // R3BEventHeader for trigger information, if needed!
    fHeader = (R3BEventHeader*)fMan->GetObject("R3BEventHeader");

    //Container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    //to be set as dependency in CMakelists.txt (in this case in the psp dir.)
    fCalPar = (R3BPspxCalPar*)FairRuntimeDb::instance()->getContainer("PspxCalPar");
    if (!fCalPar) {
        LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). Could get handle of PspxCalPar ";
        return kFATAL;
    }

    fCalPar->setChanged();

    if (!NumDet)
    {
        LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No Detectors detected."
        return kFATAL;
    }

    if (!NumStrips)
    {
        LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No Strips found."
        return kFATAL;
    }

    fMan->Register("PspxCal", "Pspx", fCalItems, kTRUE);
    return kSUCCESS;
}

void R3BPspxMapped2Cal::SetParContainers()
{
    /*** Initialize/Reads parameter file for conversion.*/
    fCalPar = (R3BPspxCalPar*)FairRuntimeDb::instance()->getContainer("PspxCalPar");
    if (!fCalPar)
    {
        LOG(ERROR) << "Could not get access to PspxCalPar-Container.";
        return;
    }
}

InitStatus R3BPspxMapped2CalPar::ReInit()
{
    /*** Initialize/Reads parameter file for conversion.*/
    SetParContainers();
    return kSUCCESS;
}

void R3BPspxMapped2CalPar::Exec(Option_t* option)
{
    /**
     * Does the conversion from Mapped to Cal level. It is called for every event.
     * Energies, which are below a channel specific threshold, will be ignored.
     * Applies (strip specific) gains to the energy entries of side 2 of the strip. This is
     * the first calibration step for the position reconstruction.
     */
    Int_t nHits->fMapped->GetEntries();

    //loop over mapped hits
    for (Int_t ii = 0; ii < nHits; ii++)
    {
        auto mapped = (R3BPspxMappedData const*)fMapped->At(ii);
        if (mapped->GetDetector() > fNumDet)
        {
            LOG(ERROR) << "R3BPspxMapped2CalPar::Exec() : more det than expected! Det: " << mapped->GetDetectorId()
            << "allowed are 1.." << fNumDet;
            continue;
        }
        if (mapped->GetNumStrips() > fNofStrips)
        {       LOG(ERROR) << "R3BPspxMapped2CalPar::Exec() : more det than expected! Det: " << mapped->GetStripId()
            << "allowed are 1.." << fNumDet;
            continue;
        }

        if (mapped->GetStrip1() == mapped->GetStrip2())
        {
            hEnergy1[mapped->GetStrip1()]->Fill(mapped->GetEnergy1());
            hEnergy2[mapped->GetStrip2()]->Fill(mapped->GetEnergy2());

            //TFit *gaus1 = new TFit();
            //if ()
            //TFit *gaus2 = new TFit();
            //if (TMath::Abs(mappedData->GetEnergy1())  &&
            //    TMath::Abs(mappedData->GetEnergy2()))
            //{ // strips
            //    Float_t energy1 = mappedData->GetEnergy1();
            //    Float_t energy2 = mappedData->GetEnergy2() * gain[d][strip - 1];
            //    Float_t time1   = mappedData->GetTime1();
            //    Float_t time2   = mappedData->GetTime2();
            //    Float_t trigger = mappedData->GetTrigger();
            //    new ((*fCalItems[d])[fCalItems[d]->GetEntriesFast()])
            //        R3BPspxCalData(strip, energy1, energy2, time1, time2, trigger);
            //}
        }
    }
}

void R3BPspxMapped2Cal::FinishEvent()
{
    fCalItems->Clear();
}

void R3BPspxMapped2Cal::FinishTask() {fCalpar->printParams();}

void R3BPspxMapped2Cal::SetUpdateRate(Int_t rate) {fUpdateRate = rate;}
void R3BPspxMapped2Cal::SetMinStats(Int_t minStats) {fMinStats = minStats;}
void R3BPspxMapped2Cal::SetNumDet(Int_t det) {fNumDet = det;}
void R3BPspxMapped2Cal::SetNumStrips(Int_t strip) {fNumStrips = strip;}

ClassImp(R3BPspxMapped2CalPar, 1)
