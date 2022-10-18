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

#ifndef R3BPSPXMAPPEDDATA_H
#define R3BPSPXMAPPEDDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Mapped level.
 * @author Ralf Plag, Ina Syndikus
 * @since January 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxMappedData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxMappedData();

    /** Standard Constructor **/
    R3BPspxMappedData(Int_t detector, Int_t face, Int_t side, Int_t strip, Int_t energy, Int_t time, Int_t trigger);

    /** Destructor **/
    virtual ~R3BPspxMappedData() {}

    void SetValue(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);
    // Getters
    inline const Int_t& GetDetector() const { return fDetector; }
    inline const Int_t& GetFace() const { return fFace; }
    inline const Int_t& GetSide() const { return fSide; }
    inline const Int_t& GetStrip() const { return fStrip; }
    inline const Int_t& GetEnergy() const { return fEnergy; }
    inline const Int_t& GetTime() const { return fTime; }
    inline const Int_t& GetTrigger() const { return fTrigger; }

  private:
    Int_t fDetector;  // Detector number, for s473 should be 3 
    Int_t fFace;  // Face number, counting from 1, Face = 1 is x and Face = 2 is y
    Int_t fSide;  // Channel number, counting from 1, one entry for each side of each strip
    Int_t fStrip;  // Channel number, counting from 1, one entry for each side of each strip
    Int_t fEnergy; // Energy/Collected charge, one entry for each side of each strip
    Int_t fTime; // Time for each hit, one entry for each side of each strip
    Int_t fTrigger; // Time for each hit, one entry for each side of each strip

  public:
    ClassDef(R3BPspxMappedData, 6)
};

#endif
