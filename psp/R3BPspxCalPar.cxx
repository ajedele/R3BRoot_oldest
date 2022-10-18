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
// -----------------------------------------------------------------
// -----           R3BPspxMappedPar header file                -----
// -----           Created 16/05/12  by I.Syndikus             -----
// -----           Modified Dec 2019 by M. Holl                -----
// -----------------------------------------------------------------

#include "R3BPspxCalPar.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParIo.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"
#include "TArrayF.h"

// -------- Standard Constructor ----------------------------------------------
R3BPspxCalPar::R3BPspxCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDet(-1)
    , fNumStrips(-1)
{
    fCalPosPar = new TArrayF(fNumDet * 2 * fNumStrips * fNumStrips); // number of detectors * number of faces * number of strips * *number of strips on other side  - pos cal par
}

// -------- Destructor --------------------------------------------------------
R3BPspxCalPar::~R3BPspxCalPar() 
{
    	clear();
	if (fCalPosPar) delete fCalPosPar;
}

// -------- Method clear ------------------------------------------------------
void R3BPspxCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// -------- Method putParams --------------------------------------------------
void R3BPspxCalPar::putParams(FairParamList* list)
{
    LOG(INFO) << "I am in R3BPspxCalPar::putParams ";
    if (!list)
    { 	
	LOG(ERROR) << "Could not find FairParamList)";    
	    return;
    }

    list->add("R3BPspxCalDetectors", fNumDet);
    list->add("R3BPspxCalStrips", fNumStrips);

    Int_t array_size = fNumDet * 2 * fNumStrips * fNumStrips;
    LOG(INFO) << "Array size: " << array_size;

    fCalPosPar->Set(array_size);
    list->add("R3BPspxPosPar", *fCalPosPar);
}

// -------- Method getParams --------------------------------------------------
Bool_t R3BPspxCalPar::getParams(FairParamList* list)
{
    LOG(INFO) << "I am in R3BPspxCalPar::getParams ";

    if (!list)
    { 	
	LOG(ERROR) << "Could not find FairParamList)";    
	    return kFALSE;
    }

    if (!list->fill("R3BPspxCalDetectors", &fNumDet)) {return kFALSE;}
    if (!list->fill("R3BPspxCalStrips", &fNumStrips)) {return kFALSE;}

    // count all entries: lines with strip info (2 entries) + lines with detector info (3 entries)
    Int_t array_size = fNumDet * 2 * fNumStrips * fNumStrips;
    LOG(INFO) << "R3BPspxCalPosPar Array Size: " << array_size;

    fCalPosPar->Set(array_size);
    if (!(list->fill("R3BPspxCalPosPar", fCalPosPar)))
    {
        LOG(WARNING) << "Could not initialize R3BPspxCalPosPar";
        return kFALSE;
    }

    return kTRUE;
}

// -------- Print Params ------------------------------------------------------
void R3BPspxCalPar::printParams()
{

    LOG(INFO) << "R3BPspxCalPosPar::printParams";
    LOG(INFO) << "fNumDetectors: " << fNumDet;
    LOG(INFO) << "fNumStrips: " << fNumStrips;

    Int_t array_size = fNumDet * 2 * fNumStrips * fNumStrips;
    LOG(INFO) << "\n***********************\nArray Size: " << array_size << "\n***********************\n";
    LOG(INFO) << "fCalPosPar: ";
    //for (Int_t i = 0; i < array_size; i++)
    //{
    //    LOG(INFO) << i << " :" << fCalPosPar->GetAt(i);
    //}
}

ClassImp(R3BPspxCalPar)
