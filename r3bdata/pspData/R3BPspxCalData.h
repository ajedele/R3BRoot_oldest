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

#ifndef R3BPSPXCALDATA_H
#define R3BPSPXCALDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Cal level.
 * @author Ralf Plag, Ina Syndikus
 * @since January 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxCalData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxCalData();

    /** Standard Constructor **/
    R3BPspxCalData(Int_t detector, Int_t stripfront, Int_t stripback, Float_t energyfront1, Float_t energyfront2, Float_t energyback1, Float_t energyback2, Float_t timefront1, Float_t timefront2, Float_t timeback1, Float_t timeback2, Int_t multfront, Int_t multback);

    /** Destructor **/
    virtual ~R3BPspxCalData() {}


    // Getters
    inline const Int_t& GetDetector() const { return fDetector; }
    inline const Int_t& GetStripFront() const { return fStripFront; }
    inline const Int_t& GetStripBack() const { return fStripBack; }
    inline const Float_t& GetEnergyFront1() const { return fEnergyFront1; }
    inline const Float_t& GetEnergyFront2() const { return fEnergyFront2; }
    inline const Float_t& GetEnergyBack1() const { return fEnergyBack1; }
    inline const Float_t& GetEnergyBack2() const { return fEnergyBack2; }
    inline const Float_t& GetTimeFront1() const { return fTimeFront1; }
    inline const Float_t& GetTimeFront2() const { return fTimeFront2; }
    inline const Float_t& GetTimeBack1() const { return fTimeBack1; }
    inline const Float_t& GetTimeBack2() const { return fTimeBack2; }
    inline const Int_t& GetMultFront() const { return fMultFront; }
    inline const Int_t& GetMultBack() const { return fMultBack; }

  private:
    Int_t fDetector; // Detector number, counting from 1
    Int_t fStripFront;    // Strip number, counting from 1
    Int_t fStripBack;    // Strip number, counting from 1
    Float_t fEnergyFront1; // Energy from both sides of hit strip after corrected for the position.
    Float_t fEnergyFront2; // Energy from both sides of hit strip after corrected for the position.
    Float_t fEnergyBack1; // Energy from both sides of hit strip after corrected for the position.
    Float_t fEnergyBack2; // Energy from both sides of hit strip after corrected for the position.
    Float_t fTimeFront1;   // Hit time from both sides of hit strip. This value is corrected relative to the trigger time
    Float_t fTimeFront2;   // Hit time from both sides of hit strip. This value is corrected relative to the trigger time
    Float_t fTimeBack1;   // Hit time from both sides of hit strip. This value is corrected relative to the trigger time
    Float_t fTimeBack2;   // Hit time from both sides of hit strip. This value is corrected relative to the trigger time
    //For Time calibration, dependent on the FW used. For Nik's FW, time has to be recentered. For CALIFA FW, time subtracted relative to the trigger.
    Int_t fMultFront;  //Multiplicity of the front side
    Int_t fMultBack;   //Multiplicity of the back side


  public:
    ClassDef(R3BPspxCalData, 9)
};

#endif
