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

#ifndef R3BPSPXCALPAR_H
#define R3BPSPXCALPAR_H

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"
#include <TObjString.h>

class FairParIo;
class FairParamList;

/**
 * Class for Parameters for Precal2Cal Conversion for PSPX detector data.
 * @author Ina Syndikus
 * @since May 12, 2016
 * Modified Dec 2019 by M.Holl
 * Modified July 2022 by A. Jedele
 */

class R3BPspxCalPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BPspxCalPar(const char* name = "PspxCalPar",
                  const char* title = "Pspx CAL parameters",
                  const char* context = "Default");

    // Destructor
    virtual ~R3BPspxCalPar();

    // Reset all parameters
    virtual void clear();

    //Method to retrieve all parameters to the standard output
    Bool_t getParams(FairParamList*);

    //Method to store all parameters using FairRunTimeDB
    void putParams(FairParamList*);

    // Print parameters
    virtual void printParams();

    // Getter & Setter
    const Int_t GetNumDet() { return fNumDet; }
    const Int_t GetNumStrips() { return fNumStrips; }
    TArrayF* GetCalPosPar() { return fCalPosPar; }

    void SetNumDet(Int_t numdet) { fNumDet = numdet; }
    void SetNumStrips(Int_t numstrips) { fNumStrips = numstrips; }
    //void SetCalPosPar(Int_t det, Int_t face, Int_t strip, Int_t backstrip, Double_t ratio) { fCalPosPar->AddAt(det,face,strip,backstrip,ratio); }
    //void SetCalPosPar(Int_t det, Int_t face, Int_t strip, Int_t backstrip, Double_t ratio) { fCalPosPar->AddAt(ratio, Index); }
    void SetCalPosPar(Int_t det, Int_t face, Int_t strip, Int_t backstrip, Double_t ratio) { fCalPosPar->AddAt(ratio, fNumStrips*fNumStrips*2*det + fNumStrips*fNumStrips*face + fNumStrips*strip + backstrip); }

  private:
    Int_t fNumDet; // number of detectors
    Int_t fNumStrips;    // number of strips per detector
    TArrayF* fCalPosPar;     // calibration parameters for each strip

    //Int_t Index = fNumStrips*fNumStrips*2*det + fNumStrips*fNumStrips*face + fNumStrips*strip + backstrip;

    R3BPspxCalPar(const R3BPspxCalPar&); //a copy constructor
    const R3BPspxCalPar& operator=(const R3BPspxCalPar&); // an  assignment operator

    ClassDef(R3BPspxCalPar, 5);
};

#endif
