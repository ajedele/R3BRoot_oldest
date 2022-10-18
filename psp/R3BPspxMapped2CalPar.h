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
// -----                   R3BPspxMapped2CalPar               -----
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl	        -----
// ----------------------------------------------------------------

#ifndef R3BPSPXMAPPED2CALPAR_H
#define R3BPSPXMAPPED2CALPAR_H

#include "FairTask.h"
#include "R3BPspxMappedData.h"
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

class R3BPspxMapped2CalPar : public FairTask
{
  public:
    /** Default Constructor **/
    R3BPspxMapped2CalPar();

    /** Standard Constructor **/
    R3BPspxMapped2CalPar(const char* name, Int_t iVerbose);

    /** Destructor **/
    virtual ~R3BPspxMapped2CalPar();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t*);

    /** Virtual method FinishTask **/
    virtual void Reset();

    /** Virtual method GetPosParamters **/
    virtual void GetPosParameters();

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();


//    /** Method SetParContainers **/
//    void SetParContainers();

    //Method to specify experiment. Very important due to different FWs used
    void SetNumExpt(Int_t);

    /** Method for setting number of detectors. 
     * No need to specify number of sides or channels 
     * - that's fixed for the PSPs **/
    void SetNumDet(Int_t);
    void SetNumStrips(Int_t);

    void SetParOutName(TString);

  private:
    R3BEventHeader* fHeader;                 // do we need that?
    TClonesArray* fMapped; /**< Arrays holding input (Mapped) data */

    UInt_t fNumDet; /**Number of detectors */
    UInt_t fNumStrips; /**Number of detectors */

    UInt_t fNumExpt; /**Number of experiment */

    TString fParOutName;

    Int_t counter_tot = 0;
    Int_t counter_tot_tot = 0;
    Int_t counter_matched_bad = 0;
    Int_t counter_timing_bad = 0;
    Int_t counter_events = 0;
    Int_t counter_interstrip = 0;
    Int_t counter_tot_interstrip = 0;
    Int_t counter_bad_interstrip = 0;
    Int_t counter_mult_mismatch = 0;
    Int_t counter_partner = 0;
    Int_t counter_good_partner = 0;
    Int_t counter_bad_time_partner = 0;
    Int_t counter_bad_energy_partner = 0;
    Int_t counter_bad_energy = 0;

    Int_t timing_index=0;
    Int_t posx_index=0;
    Int_t posy_index=0;

    Int_t maxenergy[6] = {450000, 30000, 600000, 40000, 600000, 40000};
    Double_t GainFactors[6][2][32][32] = {{{{0}}}};

    Bool_t ZombieAreAlive = kFALSE;

    Float_t GetTimeCorr_s473(Int_t time);
    Float_t GetTimeCorr_s515(Int_t time, Int_t trigger);
    Bool_t MatchedEvents(Int_t index);
    Int_t MultFace(Int_t det, Int_t face);
    Int_t InterstripEvent(Int_t ii);
    Int_t FacePartner(Int_t ii);

    Int_t fArray[10][2][2][32] = {{{{-1}}}};


    TH1F* hEnergy[10][2][2][32];
    TH1F* hTime[10][2][2][32];
    TH1F* hEnergyTot[10][2][32];
    TH1F* hInterstrip[10][2][32];
    TH1F* hInterstripTot[10][2][32];
    TH2F* hStrip1vsStrip2[10][2][32];

    TH2F* hFits[10];
    TH2F* hMult[10];
    TH1F* hTimeDiffPartner[10];
    TH2F* hHitMap[10];
    TH2F* hPosX_det[10];
    TH2F* hPosY_det[10];
    TH2F* hPosX[10][32];
    TH2F* hPosY[10][32];
    TH2F* hPosXvsStrip_det[10];
    TH2F* hPosYvsStrip_det[10];

    TH1F* hPosX_strip[10][100];
    TH1F* hPosY_strip[10][100];
       
    TDirectory *savDir = gDirectory; 
    TDirectory *Timing = savDir->mkdir("Timing"); 
    TDirectory *Energy = savDir->mkdir("Energy"); 
    TDirectory *Interstrip = savDir->mkdir("Interstrip"); 
    TDirectory *Position = savDir->mkdir("Position"); 
    TDirectory *PositionStrips = savDir->mkdir("PositionStrips"); 
    TDirectory *Fits = savDir->mkdir("Fits"); 
   
    std::ofstream outputfile;

  public:
    ClassDef(R3BPspxMapped2CalPar, 1)
};

#endif
