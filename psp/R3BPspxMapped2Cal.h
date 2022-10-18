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

#ifndef R3BPSPXMAPPED2CALTEST_H
#define R3BPSPXMAPPED2CALTEST_H

#include "FairTask.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxCalData.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"

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
 * Modified July 2022 by A. Jedele
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

    /** Virtual method FinishTask **/
    virtual void SetParContainers();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

//    /** Virtual method FinishTask **/
//    virtual void Reset();

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();


    /** Method for setting number of detectors. 
     * No need to specify number of sides or channels 
     * - that's fixed for the PSPs **/
    void SetNumExpt(Int_t);
    void SetNumDet(Int_t);
    void SetNumStrips(Int_t);

  private:
    R3BEventHeader* fHeader;                 // do we need that?
    TClonesArray* fMapped; /**< Arrays holding input (Mapped) data */
    TClonesArray* fCal; /**< Arrays holding output (Cal) data */

    R3BPspxCalPar* fCalPar; /**< Parameter instance holding thresholds and gains for position correction */

    UInt_t fNumExpt; /**Number of detectors */
    UInt_t fNumDet; /**Number of detectors */
    UInt_t fNumStrips; /**Number of detectors */

    
    Int_t maxenergy[6] = {450000, 30000, 600000, 40000, 600000, 40000};

    Bool_t ZombieAreAlive = kFALSE;

    Float_t GetTimeCorr_s473(Int_t time);
    Float_t GetTimeCorr_s515(Int_t time, Int_t trigger);
    Bool_t MatchedEvents(Int_t index);
    Int_t MultFace(Int_t det, Int_t face);
    Int_t FacePartner(Int_t ii);

    Int_t fArray[10][2][2][32] = {{{{-1}}}};
    Double_t GainFactors[6][2][32][32] = {{{{0.}}}};

    TH2F* hTimeFront[10][2][32];
    TH2F* hTimeBack[10][2][32];
    TH2F* hEnergyFront[10][2][32];
    TH2F* hEnergyCorrFront[10][2][32];
    TH2F* hEnergyBack[10][2][32];
    TH2F* hEnergyCorrBack[10][2][32];

    Int_t counter=0; Int_t counter_total=0;

  public:
    ClassDef(R3BPspxMapped2Cal, 1)
};

#endif
