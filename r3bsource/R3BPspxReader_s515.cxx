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
//#include "FairRunAna.h"
//#include "FairRuntimeDb.h"

#include "R3BPspxReader_s515.h"
#include "R3BPspxMappedData.h"
//#include "R3BPspxMappedPar.h"

#include "TClonesArray.h"
#include "ext_data_struct_info.hh"
#include <iostream>

extern "C"
{
#include "ext_data_client.h"
#include "ext/ext_h101_pspx_s515.h"
}

#define LENGTH(x) (sizeof x / sizeof *x)
#define MAX_PSPX_DETECTORS (sizeof data->PSPX / sizeof data->PSPX[0])

R3BPspxReader_s515::R3BPspxReader_s515(EXT_STR_h101_PSPX* data, UInt_t offset)
    : R3BReader("R3BPspxReader_s515")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BPspxMappedData")) // number of faces of detectors
{
}

R3BPspxReader_s515::~R3BPspxReader_s515()
{
    if (fArray) {delete fArray;}
}

/**
 * Initialize output data. Read input data.
 */
Bool_t R3BPspxReader_s515::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(INFO) << "R3BPspxReader_s515::Init";
    EXT_STR_h101_PSPX_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_PSPX, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    //Register output data array in tree 
    FairRootManager::Instance()->Register("PspxMapped", "PSPX", fArray, kTRUE);
    //Reset();


    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    // Register output array in tree
    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;
    for (Int_t dd = 0; dd < MAX_PSPX_DETECTORS; dd++)
    {
        for (Int_t ff = 0; ff < 2; ff++)
        {
		for (Int_t ss = 0; ss < 2; ss++)
		{
			data->PSPX[dd].F[ff].S[ss].E = 0;
			data->PSPX[dd].F[ff].S[ss].TM = 0;
			data->PSPX[dd].F[ff].S[ss].TL = 0;
			data->PSPX[dd].F[ff].S[ss].W = 0;
			data->PSPX[dd].F[ff].S[ss].S = 0;
              }
        }
    }
    
//    multPerStrip[MAX_PSPX_DETECTORS*4*16] = {0};

    return kTRUE;
}

/**
 * Does the unpacking to Mapped level. It is called for every event.
 * Converts plain raw data to multi-dimensional array.
 * Ignores energies with an error message.
 */
Bool_t R3BPspxReader_s515::Read()
{
    EXT_STR_h101_PSPX_onion* data = (EXT_STR_h101_PSPX_onion*)fData;

    // loop over all detectors
    for (int32_t det = 0; det < MAX_PSPX_DETECTORS; det++)
    {
        // loop over faces
        for (int32_t face = 0; face < 2; face++)
        {
            std::vector<R3BPspxMappedData*> datas(LENGTH(data->PSPX[0].F[0].S[0].EI));
            // loop over strip sides
            for (int32_t side = 0; side < 2; side++)
            {
                auto const& dfs = data->PSPX[det].F[face].S[side];
                uint32_t numHitsE = dfs.E;
                // loop over hits
                for (uint32_t i = 0; i < numHitsE; i++)
                {
      		    int32_t strip        = dfs.EI[i]; // counting from 1 to max number of channels for an detector
                    int32_t energy       = dfs.Ev[i];
                    int32_t timestrip    = dfs.TLI[i];
                    int32_t time         = dfs.TLv[i];
                    int32_t triggerstrip = dfs.SI[i];
                    int32_t trigger      = dfs.Sv[i];
                    if (energy < 0)
                        energy = -1. * energy; // make sure energy values are positive. Necessary for compatibilty with
                                               // GSI Febex firmware

                    new ((*fArray)[fArray->GetEntriesFast()])
                    R3BPspxMappedData(det+1, face+1, side+1, strip, energy, time, trigger);
                    //datas[strip - 1]->SetValue(det, side, strip, energy, time, trigger);
		}
            }
        }
    }
    return kTRUE;
}

void R3BPspxReader_s515::Reset()
{
    fArray->Clear();
}

ClassImp(R3BPspxReader_s515)
