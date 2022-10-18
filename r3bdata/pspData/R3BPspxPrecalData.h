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

#ifndef R3BPSPXPRECALDATA_H
#define R3BPSPXPRECALDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Precal level.
 * @author Ina Syndikus
 * @since March 2017
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxPrecalData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxPrecalData();

    /** Standard Constructor **/
    R3BPspxPrecalData(Int_t detector, Int_t strip, Float_t energy1, Float_t energy2, Float_t time1, Float_t time2);

    /** Destructor **/
    virtual ~R3BPspxPrecalData() {}

    // Getters
    inline const Int_t& GetDetector() const { return fDetector; }
    inline const Int_t& GetStrip() const { return fStrip; }
    inline const Float_t& GetEnergy1() const { return fEnergy1; }
    inline const Float_t& GetEnergy2() const { return fEnergy2; }
    inline const Float_t& GetTime1() const { return fTime1; }
    inline const Float_t& GetTime2() const { return fTime2; }

  private:
    Int_t fDetector;  // Detector number, counting from 1 */
    Int_t fStrip;     // Strip number, counting from 1 */
    Float_t fEnergy1; // Energy/Collected charge for side 1 of the strip. This value needs no gain correction
                      // (Mapped2Precal Position Calibration).
    Float_t fEnergy2; // Energy/Collected charge for side 2 of  the strip. This value is gain corrected (Mapped2Precal
                      // Position Calibration).
    Float_t fTime1;   // Hit time for side 1 of the strip. This value is corrected relative to the internal master trigger time.
    Float_t fTime2;   // Hit time for side 2 of  the strip. This value is corrected relative to the internal master trigger time.

  public:
    ClassDef(R3BPspxPrecalData, 4)
};

#endif
