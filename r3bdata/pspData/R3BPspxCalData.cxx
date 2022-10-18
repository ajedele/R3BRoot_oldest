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

#include "R3BPspxCalData.h"

R3BPspxCalData::R3BPspxCalData()
    : fDetector(-1)
    , fStripFront(-1)
    , fStripBack(-1)
    , fEnergyFront1(-1)
    , fEnergyFront2(-1)
    , fEnergyBack1(-1)
    , fEnergyBack2(-1)
    , fTimeFront1(-1)
    , fTimeFront2(-1)
    , fTimeBack1(-1)
    , fTimeBack2(-1)
    , fMultFront(-1)
    , fMultBack(-1)
{
}

R3BPspxCalData::R3BPspxCalData(Int_t detector, Int_t stripfront, Int_t stripback, Float_t energyfront1, Float_t energyfront2, Float_t energyback1, Float_t energyback2, Float_t timefront1, Float_t timefront2, Float_t timeback1, Float_t timeback2, Int_t multfront, Int_t multback)
    : fDetector(detector)
    , fStripFront(stripfront)
    , fStripBack(stripback)
    , fEnergyFront1(energyfront1)
    , fEnergyFront2(energyfront2)
    , fEnergyBack1(energyback1)
    , fEnergyBack2(energyback2)
    , fTimeFront1(timefront1)
    , fTimeFront2(timefront2)
    , fTimeBack1(timeback1)
    , fTimeBack2(timeback2)
    , fMultFront(multfront)
    , fMultBack(multback)
{
}

ClassImp(R3BPspxCalData)
