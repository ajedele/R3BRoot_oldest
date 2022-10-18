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

#ifndef R3BPSPXDIGITIZER_H
#define R3BPSPXDIGITISER_H 1

#include "FairTask.h"
#include "R3BPspxDigi.h"
#include "R3BPspxDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BPspxDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BPspxDigitizer();

    /** Destructor **/
    ~R3BPspxDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BPspxDigi* AddHit(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e);

  protected:
    TClonesArray* fPspxPoints;
    TClonesArray* fPspxMCTrack;
    TClonesArray* fPspxDigi;

    // Parameter class
    R3BPspxDigiPar* fPspxDigiPar;

    //- Control Hitograms

    //   TH1F* PSPXhis;

    Int_t eventNoPspx;

  private:
    virtual void SetParContainers();

    ClassDef(R3BPspxDigitizer, 1);
};

#endif
