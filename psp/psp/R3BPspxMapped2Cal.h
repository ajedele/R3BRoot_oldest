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

// ----------------------------------------------------------------
// -----                   R3BPspxMapped2Cal               -----
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl	        -----
// ----------------------------------------------------------------

#ifndef R3BPSPXMAPPED2PRECAL_H
#define R3BPSPXMAPPED2PRECAL_H

#include "FairTask.h"

class TClonesArray;
class R3BEventHeader;
class R3BPspxCalPar;

/**
 * Class to convert Mapped data to Cal data for PSPX detector data.
 * Thresholds are applied to signal from both sides of each strip
 * Signal from side 2 of each strip is multiplied by gain for position calibration
 * @author Ina Syndikus
 * @since March 13, 2016
 * Modified Dec 2019 by M.Holl
 * Modified April 2021 by J.L.Rodriguez
 */

class R3BPspxMapped2Cal : public FairTask
{
  public:
    /** Default Constructor **/
    R3BPspxMapped2Cal();

    /** Standard Constructor **/
    R3BPspxMapped2Cal(const char* name, Int_t iVerbose);

    /** Destructor **/
    virtual ~R3BPspxMapped2Cal();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Method SetParContainers **/
    void SetParContainers();

    /** Method for setting number of detectors. 
     * No need to specify number of sides or channels 
     * - that's fixed for the PSPs **/
    void SetNumDet(Int_t);
    void SetNumStrips(Int_t);
    void SetMinStats(Int_t);
    void SetUpdateRate(Int_t);

  private:
    void SetParameters();

    R3BEventHeader* fHeader;                 // do we need that?
    std::vector<TClonesArray*> fMapped; /**< Arrays holding input (Mapped) data */
    std::vector<TClonesArray*> fCal; /**< Arrays holding output (Cal) data */

    R3BPspxCalPar* fCalPar; /**< Parameter instance holding thresholds and gains for position correction */
    std::vector<std::vector<Float_t>> gain;
    std::vector<std::vector<Int_t>> threshold1;
    std::vector<std::vector<Int_t>> threshold2;

    UInt_t fNumDet; /**Number of detectors */
    UInt_t fNumStrips; /**Number of detectors */
    Int_t fUpdateRate; /**Rate maximum */
    Int_t fMinStats; /**Minimum Statistics */
    Int_t fTpat; /**Minimum Statistics */

  public:
    ClassDef(R3BPspxMapped2Cal, 3)
};

#endif
