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

#include "R3BPspxMappedData.h"

R3BPspxMappedData::R3BPspxMappedData()
    : fDetector{ -1 }
    , fFace{ -1 }
    , fSide{ -1 }
    , fStrip{ -1 }
    , fEnergy{-1 }
    , fTime{ -1 }
    , fTrigger{ -1 }
{
}

R3BPspxMappedData::R3BPspxMappedData(Int_t detector, Int_t face, Int_t side, Int_t strip, Int_t energy, Int_t time, Int_t trigger)
    : fDetector{ detector }
    , fFace{ face }
    , fSide{ side }
    , fStrip{ strip }
    , fEnergy{ energy }
    , fTime{ time }
    , fTrigger{ trigger }
{
}

void R3BPspxMappedData::SetValue(Int_t detector, Int_t face, Int_t side, Int_t strip, Int_t energy, Int_t time, Int_t trigger)
{
    fDetector = detector;
    fFace = face;
    fSide = side;
    fStrip = strip;
    fEnergy = energy;
    fTime = time;
    fTrigger = trigger;
}

ClassImp(R3BPspxMappedData)
