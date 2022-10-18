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
// -----                      R3BPspxPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BPspxDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BPspxDigi::R3BPspxDigi()
    : Ps03mul(0)
{
}

R3BPspxDigi::R3BPspxDigi(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e)
    : R3BHit(0, psp3x, psp3y, psp3e, 0.)
    , Ps03mul(psp3mul)
{
}

// -----   Destructor   ----------------------------------------------------
R3BPspxDigi::~R3BPspxDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BPspxDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BPspxDigi)
