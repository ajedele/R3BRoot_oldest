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
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl          -----
// -----              Modified  Jul 2022  by A. Jedele        -----
// ----------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxMapped2Cal.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxCalData.h"
#include "R3BPspxCalPar.h"

#include "TClonesArray.h"
#include <iostream>
#include <limits>

R3BPspxMapped2Cal::R3BPspxMapped2Cal()
    : FairTask("PspxMapped2Cal", 1)
    , fNumDet(0)
//    , fMappedItems()
//    , fCalItems()
{
}

R3BPspxMapped2Cal::R3BPspxMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumDet(0)
    , fNumExpt(0)
    , fNumStrips(0)
    , fMapped()
//    , fMappedItems()
//    , fCalItems()
{
}

R3BPspxMapped2Cal::~R3BPspxMapped2Cal()
{

//    LOG(INFO) << "R3BPspxMapped2Cal: Delete instance";
//    for (Int_t i = 0; i < fMappedItems.size(); i++)
//    {
//        delete fMappedItems[i];
//    }
//    for (Int_t i = 0; i < fCalItems.size(); i++)
//    {
//        delete fCalItems[i];
//    }
}

void R3BPspxMapped2Cal::SetParameters()
{
    LOG(INFO) << "R3BPspxMapped2Cal::SetParameters()";
    //--- Parameter Container ---
    //Int_t nDet = fCalPar->GetNumDetectors(); // Number of Detectors/Faces
    Int_t nDet = fNumDet; // Number of Detectors/Faces
    printf("Number of Detectors: %d\n",nDet);
    LOG(INFO) << nDet;
    gain.resize(nDet);
    threshold1.resize(nDet);
    threshold2.resize(nDet);
    for (Int_t d = 0; d < nDet; d++)
    {
        Int_t nStrips = fNumStrips; // Number of Strips
        gain[d].resize(nStrips);
        threshold1[d].resize(nStrips);
        threshold2[d].resize(nStrips);
        Int_t parOffset =
            d * nStrips * 4 +
            (d + 1) * 3; // Position in parameter list. 4 parameters per strip + 3 "header" parameters per detector.
        for (Int_t s = 0; s < nStrips; s++)
        {
            //TArrayF par = fCalPar->GetCalPar(); // Array with the parameters
            //gain[d][s] = par.At(parOffset + 1);
            //threshold1[d][s] = par.At(parOffset + 2);
            //threshold2[d][s] = par.At(parOffset + 3);
            //parOffset += 4; // move to next line in parameter file.
            //LOG(INFO) << "Det: " << d << "\tstr: " << s << "\tgain: " << gain[d][s] << "\tthr1: " << threshold1[d][s]
            //          << "\tthr2: " << threshold2[d][s];
        }
    }
}

InitStatus R3BPspxMapped2Cal::Init()
{
    /**
     * Initialize output data. Read input data and parameters.
     * The parameters get saved in dedicated arrays.
     * Print parameters, if verbosity is set to INFO.
     */

    LOG(INFO) << "R3BPspxMapped2Cal::Init()";
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(ERROR) << "R3BPspxMapped2Cal::Init() Root-manager not found.";
        return kFATAL;
    }

    // R3BEventHeader for trigger information, if needed!
    fHeader = (R3BEventHeader*)fMan->GetObject("R3BEventHeader");

    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    // Figure out how many detectors were registered by the reader
    for (Int_t d = 0; d < fCalPar->GetNumDetectors(); d++)
    {
        TClonesArray* tmp[2];

        for (Int_t f = 0; f < 2; f++)
        {
            tmp[f] = (TClonesArray*)fMan->GetObject(Form("Pspx%d_%cMapped", d + 1, xy[f])); // = branch name in TTree
        }
        if (tmp[0] == NULL && tmp[1] == NULL)
        {
            if (d == 0)
            {
                LOG(ERROR) << "R3BPspxMapped2Cal::Init() Couldn't get handle on PSPX-mapped items.";
                return kFATAL;
            }
            break;
        }
        for (Int_t f = 0; f < 2; f++)
        {
//            fMappedItems.push_back(tmp[f]);
//            fCalItems.push_back(new TClonesArray("R3BPspxCalData"));
//            fMan->Register(
//                Form("Pspx%d_%cCal", d + 1, xy[f]), Form("Pspx%d_%c", d + 1, xy[f]), fCalItems.back(), kTRUE);
        }
    }

//    SetParameters();
    return kSUCCESS;
}

void R3BPspxMapped2Cal::SetParContainers()
{
    /*** Initialize/Reads parameter file for conversion.*/

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
//        LOG(ERROR) << "R3BPspxMapped2Cal::FairRuntimeDb not opened!";
//        return;
    }
    else
    {
//        LOG(INFO) << "R3BPspxMapped2Cal::SetParContainers()";
    }

    fCalPar = (R3BPspxCalPar*)rtdb->getContainer("R3BPspxCalPar");
    if (!fCalPar)
    {
//        LOG(ERROR) << "R3BPspxMapped2Cal::Could not get access to R3BPspxCalPar-Container.";
//        return;
    }

    fCalPar->printParams();
}

InitStatus R3BPspxMapped2Cal::ReInit()
{
    /*** Initialize/Reads parameter file for conversion.*/

    LOG(INFO) << "R3BPspxMapped2Cal::ReInit()";
    /*
        fCalPar = (R3BPspxCalPar*)FairRuntimeDb::instance()->getContainer("R3BPspxCalPar");

        if (!fCalPar)
        {
            LOG(ERROR) << "Could not get access to R3BPspxCalPar-Container.";
            return kFATAL;
        }*/

    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

void R3BPspxMapped2Cal::Exec(Option_t* option)
{
    /**
     * Does the conversion from Mapped to Cal level. It is called for every event.
     * Energies, which are below a channel specific threshold, will be ignored.
     * Applies (strip specific) gains to the energy entries of side 2 of the strip. This is
     * the first calibration step for the position reconstruction.
     */
/*
    for (Int_t d = 0; d < fMappedItems.size(); d++)
    {
        if (!fMappedItems[d])
        {
            printf("Cannot access PSPX%d_%d mapped items\n", (d / 2) + 1, (d % 2) + 1);
            return;
        }
        if (fCalPar->GetNumStrips().At(d) == 0)
            continue;
        Int_t nMapped = fMappedItems[d]->GetEntries();

        for (Int_t i = 0; i < nMapped; i++)
        {

            R3BPspxMappedData* mappedData = (R3BPspxMappedData*)fMappedItems[d]->At(i);

            // get rid of error message from Febex (GSI firmware)
            //if (mappedData->GetEnergy1() == 3075811 || mappedData->GetEnergy1() == 3075810)
            //    continue;
            //if (mappedData->GetEnergy2() == 3075811 || mappedData->GetEnergy2() == 3075810)
            //    continue;
            
            if (mappedData->GetStrip1() == mappedData->GetStrip2())
            {
                Int_t strip = mappedData->GetStrip1();
                if (TMath::Abs(mappedData->GetEnergy1()) > threshold1[d][strip - 1] &&
                    TMath::Abs(mappedData->GetEnergy2()) > threshold2[d][strip - 1])
                { // strips
                    Float_t energy1 = mappedData->GetEnergy1();
                    Float_t energy2 = mappedData->GetEnergy2() * gain[d][strip - 1];
                    Float_t time1   = mappedData->GetTime1();
                    Float_t time2   = mappedData->GetTime2();
                    Float_t trigger = mappedData->GetTrigger();
                    new ((*fCalItems[d])[fCalItems[d]->GetEntriesFast()])
                        R3BPspxCalData(strip, energy1, energy2, time1, time2, trigger);
                }
            }
        }
    }
*/}

void R3BPspxMapped2Cal::FinishEvent() {}

void R3BPspxMapped2Cal::FinishTask()
{
//    for (Int_t i = 0; i < fMappedItems.size(); i++)
//    {
//        fMappedItems[i]->Clear();
//    }
//    for (Int_t i = 0; i < fCalItems.size(); i++)
//    {
//        fCalItems[i]->Clear();
//    }
}


void R3BPspxMapped2Cal::SetNumDet(Int_t det)
{
    fNumDet = det;
}

ClassImp(R3BPspxMapped2Cal)
