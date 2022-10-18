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
// -----                   R3BPspxCal2HitPar               -----
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl	        -----
// ----------------------------------------------------------------

#ifndef R3BPSPXCAL2HITPAR_H
#define R3BPSPXCAL2HITPAR_H

#include "FairTask.h"
#include "R3BPspxMappedData.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TCanvas.h"

class TClonesArray;
class R3BEventHeader;
class R3BPspxHitPar;

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

class R3BPspxCal2HitPar : public FairTask
{
  public:
    /** Default Constructor **/
    R3BPspxCal2HitPar();

    /** Standard Constructor **/
    R3BPspxCal2HitPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BPspxCal2HitPar();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t*);

    /** Virtual method FinishTask **/
    virtual void Reset();

    /** Virtual method FinishTask **/
    virtual void FinishTask();


    /** Virtual method GetEnergyCorrection **/
    void GetEnergyCorrection();

    /** Method for setting number of detectors. 
     * No need to specify number of sides or channels 
     * - that's fixed for the PSPs **/
    void SetNumDet(Int_t);
    void SetNumStrips(Int_t);

    void SetParOutName(TString);


  private:
    TClonesArray* fCal;    

    Int_t counter_tot = 0;
    Int_t counter_events = 0;

    UInt_t fNumDet; /** Number of detectors */
    UInt_t fNumStrips; /** Number of detectors */

    TString fParOutName;

    TH1F* RawEnergyFront[6][32][32];
    TH1F* RawEnergyBack[6][32][32];
    TH2D* RawEnergyTimeFront[6][32][32];
    TH2D* RawEnergyTimeBack[6][32][32];
    TH2F* RawEnergyFront2D[6][32];
    TH2F* RawEnergyBack2D[6][32];
      
   TCanvas *GausFitCan[6]; 
    //TDirectory *savDir = gDirectory; 

    Float_t fArray[6][32][32] = {{{0}}};  

    std::ofstream outputfile;

  public:
    ClassDef(R3BPspxCal2HitPar, 2)
};

#endif
