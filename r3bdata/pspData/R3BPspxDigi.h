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

// -------------------------------------------------------------------------
// -----                      R3BPspxDigi header file                  -----
// -----                  Created 28/06/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BPspxDigi.h
 **/

#ifndef R3BPSPXDIGI_H
#define R3BPSPXDIGI_H

#include "R3BHit.h"
#include "TVector3.h"

class R3BPspxDigi : public R3BHit
{

  public:
    /** Default constructor **/
    R3BPspxDigi();
    R3BPspxDigi(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e);

    /** Copy constructor **/
    R3BPspxDigi(const R3BPspxDigi& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BPspxDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetPspx3mul(Int_t mul) { Ps03mul = mul; }
    Double_t GetPspx3mul() { return Ps03mul; }

  protected:
    Int_t Ps03mul;

    ClassDef(R3BPspxDigi, 2)
};

#endif
