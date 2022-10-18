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

#include "FairLogger.h"

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BPspxReader.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxMappedPar.h"
#include "TClonesArray.h"
#include "ext_data_struct_info.hh"

extern "C"
{
#include "ext_data_client.h"
#include "ext/ext_h101_pspx.h"
}

#define LENGTH(x) (sizeof x / sizeof *x)

R3BPspxReader::R3BPspxReader(EXT_STR_h101_PSPX* data, UInt_t offset)
    : R3BReader("R3BPspxReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fLogger(FairLogger::GetLogger())
    //, fMappedItems(2 * LENGTH(((EXT_STR_h101_PSPX_onion*)data)->PSPX)) // number of faces of detectors
    , fMappedItems(new TClonesArray("R3BPspxMappedData")) // number of faces of detectors
{
//    EXT_STR_h101_PSPX_onion* data_o = (EXT_STR_h101_PSPX_onion*)fData;
//    for (Int_t d = 0; d < 2 * LENGTH(data_o->PSPX); d++)
//    {
//        fMappedItems[d] = new TClonesArray("R3BPspxMappedData");
//    }
//    printf("Length: %lu\n", LENGTH(data_o->PSPX));
//    LOG(INFO) << "R3BPspxReader: Created " << 2 * LENGTH(data_o->PSPX) << " detectors.";
}

R3BPspxReader::~R3BPspxReader()
{

    if (fMappedItems)
    {
	    delete fArray;
    }
//    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;
//    for (Int_t d = 0; d < 2 * LENGTH(data->PSPX); d++)
//    {
//        delete fMappedItems[d];
//    }
}

/**
 * Initialize output data. Read input data.
 */
Bool_t R3BPspxReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(INFO) << "R3BPspxReader::Init";
    EXT_STR_h101_PSPX_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_PSPX, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }
    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    // Register output array in tree
    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;
    for (Int_t dd = 0; dd < MAX_TOFD_PLANES; dd++)
    {
        for (Int_t ff = 0; ff < 2; ff++)
        {
		for (Int_t ss = 0; ss < 2; ss++)
		{
			data->PSPX[dd].F[ff].S[ss].E = 0;
			data->PSPX[dd].F[ff].S[ss].T = 0;
			data->PSPX[dd].F[ff].S[ss].TT = 0;


		//            // Register output array in tree
//            if (!fOnline)
//            {
//                FairRootManager::Instance()->Register(Form("Pspx%d_%cMapped", dd + 1, xy[ff]),
//                                                      Form("Pspx%d_%c", dd + 1, xy[ff]),
//                                                      fMappedItems[2 * dd + ff],
//                                                      kTRUE);
//            }
//            else
//            {
//                FairRootManager::Instance()->Register(Form("Pspx%d_%cMapped", dd + 1, xy[ff]),
//                                                      Form("Pspx%d_%c", dd + 1, xy[ff]),
//                                                      fMappedItems[2 * dd + ff],
//                                                      kFALSE);
//            }
//            LOG(INFO) << "Registered Pspx" << d + 1 << "_" << xy[f];
              }
        }
    }
    Reset();
    memset(fData, 0, sizeof*fData);
    return kTRUE;
}

/**
 * Does the unpacking to Mapped level. It is called for every event.
 * Converts plain raw data to multi-dimensional array.
 * Ignores energies with an error message.
 */
Bool_t R3BPspxReader::Read()
{
    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;

#if 0
    // this is the data structure we have to read:
    struct {
      struct {
        struct {
          uint32_t E;
          uint32_t EI[32 /* E */];
          uint32_t Ev[32 /* E */];
        } S[2];
      } F[2];
    } PSPX[N];

    F = Face (front/back of detector)
    S = Side (end of strip)
#endif

    // loop over all detectors
    for (Int_t det = 0; det < LENGTH(data->PSPX); det++)
    {
        // loop over faces
        for (Int_t face = 0; face < 2; face++)
        {
            std::vector<R3BPspxMappedData*> datas(LENGTH(data->PSPX[0].F[0].S[0].EI));
            // loop over strip sides
            for (Int_t side = 0; side < 2; side++)
            {
                auto const& dfs = data->PSPX[det].F[face].S[side];
                uint32_t numChannelsE = dfs.E;
                uint32_t numChannelsT = dfs.T;
                // loop over channels
                for (Int_t i = 0; i < numChannelsE; i++)
                {
                    int32_t strip   = dfs.EI[i]; // counting from 1 to max number of channels for an detector
                    int32_t energy  = dfs.Ev[i];
                    int32_t time    = dfs.Tv[i];
		    int32_t trigger;
                    if (numChannelsE < 16) {trigger = dfs.TT[0];}
                    else {trigger = dfs.TT[1];}
                   // int32_t trigger = dfs.TT;
                    if (energy < 0)
                        energy = -1. * energy; // make sure energy values are positive. Necessary for compatibilty with
                                               // GSI Febex firmware

                    // if (!datas[strip-1]) datas[strip-1] = new ((*fMappedItems[d])[fMappedItems[d]->GetEntriesFast()])
                    // R3BPspxMappedData(f+1,strip);
                    if (!datas[strip - 1])
                        datas[strip - 1] = new ((*fMappedItems[2 * det + face])[fMappedItems[2 * det + face]->GetEntriesFast()])
                            R3BPspxMappedData();
                    // assert(-1 == datas[strip-1]->GetEnergy(s));
                    datas[strip - 1]->SetValue(det, side, strip, energy, time, trigger);
                }
            }
        }
    }
    return kTRUE;
}

void R3BPspxReader::Reset()
{
    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;
    for (Int_t det = 0; det < 2 * LENGTH(data->PSPX); det++)
    {
        fMappedItems[det]->Clear();
    }
}

ClassImp(R3BPspxReader)
