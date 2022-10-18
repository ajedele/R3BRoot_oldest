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
// -----                 R3BPspxCal2Hit              -----
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl          -----
// -----              Modified  Jul 2022  by A. Jedele        -----
// ----------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxCal2Hit.h"
#include "R3BPspxCalData.h"
#include "R3BPspxCalData.h"
#include "R3BPspxHitPar.h"

#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include <limits>

//ROOT headers
#include "TH1F.h" 
#include "TH2F.h" 
#include "TMath.h" 
#include "TProfile.h" 
#include "TRandom.h" 
#include "TFile.h"
#include "TGraphErrors.h"

R3BPspxCal2Hit::R3BPspxCal2Hit()
    : FairTask("PspxCal2Hit", 1)
{
}

R3BPspxCal2Hit::R3BPspxCal2Hit(const char* name, Int_t iVerbose)
    	: FairTask(name, iVerbose)
    	, fNumExpt(0)
    	, fNumDet(0)
	, fNumStrips(0)
	, fHitPar(nullptr)
	, fCal(nullptr)
{
}

R3BPspxCal2Hit::~R3BPspxCal2Hit()
{
	delete fHitPar;
	delete fCal;
}

InitStatus R3BPspxCal2Hit::Init()
{
	//FairRootManager
	FairRootManager* fMan = FairRootManager::Instance();
	if (!fMan) { LOG(ERROR) << "R3BPspxCal2Hit::Init(). No FairRootManager found"; return kFATAL;}
	
	//Get Cal Data
	fCal = (TClonesArray*)fMan->GetObject("PspxCal");
	if (!fCal) {LOG(ERROR) << "R3BPspxCal2Hit::Init(). No PspxCalData found."; return kFATAL;}
	
	// R3BEventHeader for trigger information, if needed!
	fHeader = (R3BEventHeader*)fMan->GetObject("R3BEventHeader");
	
	//Container needs to be created in tcal/R3BCalContFact.cxx AND R3BCal needs
	//to be set as dependency in CMakelists.txt (in this case in the psp dir.)
	FairRuntimeDb* rtdb = FairRuntimeDb::instance();
	if (!rtdb) {return kFATAL;}
	
	//Get the Cal Par container
	fHitPar = (R3BPspxHitPar*)rtdb->getContainer("PspxHitPar");
	if (!fHitPar) {
	    LOG(ERROR) << "R3BPspxCal2Hit::Init(). Could get handle of PspxHitPar ";
	    return kFATAL;
	}
	
	if (!fNumDet)
	{
	    LOG(ERROR) << "R3BPspxCal2Hit::Init(). No Detectors detected.";
	    return kFATAL;
	}
	
	if (!fNumStrips)
	{
	    LOG(ERROR) << "R3BPspxCal2Hit::Init(). No Strips found.";
	    return kFATAL;
	}
	//fMan->Register("PspxHitPar", "PSPX", fCalItems, kTRUE);


	//printf("\n\nNumDet: %d, NumStrips: %d\n\n\n",fNumDet,fNumStrips);
	//TString Name;
	//for (Int_t ii = 0; ii < fNumDet; ii++)
	//{
    	//	for (Int_t jj = 0; jj < 2; jj++)
    	//	{
	//   		for  (Int_t ll = 0; ll<fNumStrips; ll++)
    	//		{
	//		}
	//	}
	//}

    return kSUCCESS;
}

InitStatus R3BPspxCal2Hit::ReInit() {return kSUCCESS;}

void R3BPspxCal2Hit::Exec(Option_t* option)
{
	Int_t nHits = fCal->GetEntries();
	if (nHits > 0)
	{
		for (Int_t ii = 0; ii < nHits; ii++)
    		{
			auto cal = (R3BPspxCalData const*)fCal->At(ii);
			if (cal == NULL) {ZombieAreAlive = kTRUE; printf("I'm am Zombie\n"); break;}  // removes bad events where -nan is filled into the fCal array
			fArray[cal->GetDetector()-1][cal->GetFace()-1][cal->GetSide()-1][cal->GetStrip()-1] = ii; 
		}

		for (Int_t ii=0; ii<nHits; ii++)
		{
			if (ZombieAreAlive == kTRUE) {break;}

			auto cal = (R3BPspxCalData const*)fCal->At(ii);
			if (MatchedEvents(ii) == kFALSE) {continue;}

			//Used to get the array index for the opposite side and opposite face of a detector
			Int_t OtherSide = 0; Int_t OtherFace = 0;
			if (cal->GetSide() == 1) {OtherSide = 2;}
			else if (cal->GetSide() == 2) {OtherSide = 1;}
			if (cal->GetFace() == 1) {OtherFace = 2;}
			else if (cal->GetFace() == 2) {OtherFace = 1;}

			timediff = cal->GetTime1() - cal->GetTime2();

			if (timediff > -10 && timediff < 10)
			{

			}
		}
	}
}

void R3BPspxCal2Hit::FinishEvent() 
{
	for (Int_t ii=0; ii < fNumDet; ii++)
	{
		for (Int_t jj=0; jj < 2; jj++)
		{
			for (Int_t kk=0; kk < 2; kk++)
			{
				for (Int_t ll=0; ll < fNumStrips; ll++)
				{
					fArray[ii][jj][kk][ll]=-1;
				}
			}
		}
	}

}

void R3BPspxCal2Hit::FinishTask()
{
	for (Int_t ii=0; ii < fNumDet; ii++)
	{
		for (Int_t jj=0; jj < 2; jj++)
		{
			for (Int_t ll=0; ll < fNumStrips; ll++)
			{
				hEnergy[ii][jj][ll]->Write();
				hEnergyCorr[ii][jj][ll]->Write();
			}
		}
	}
	fHitPar->setChanged();
}

Float_t R3BPspxCal2Hit::GetTimeCorr_s473(Int_t time)
{
	//corrects for the timing offset
	Int_t corr_time = time - 2048;
	if (corr_time < 0) {corr_time = 4096 + corr_time;}
	return corr_time;
}

Float_t R3BPspxCal2Hit::GetTimeCorr_s515(Int_t time, Int_t trigger)
{
	//corrects for the timing offset
	Int_t corr_time = time - trigger;
	return corr_time;
}

Bool_t R3BPspxCal2Hit::MatchedEvents(Int_t index) 
{
	Int_t index2 = -1, index3 = -1, index4 = -1;
	auto cal = (R3BPspxCalData const*)fCal->At(index);

        //printf("Event Number: %d, Det: %d, Face: %d, Side: %d, Strip: %d, Index: %d\n",counter_tot,cal->GetDetector(),cal->GetFace(),cal->GetSide(),cal->GetStrip(), fArray[cal->GetDetector()-1][cal->GetFace()-1][cal->GetSide()-1][cal->GetStrip()-1]);
	if (fCal->GetEntries() < 4) {return kFALSE;}
	Int_t OtherFace = -1;
	Int_t OtherSide = -1;
	if (cal->GetFace() == 1) {OtherFace = 2;}
	else if (cal->GetFace() == 2) {OtherFace = 1;}
	if (cal->GetSide() == 1) {OtherSide = 2;}
	else if (cal->GetSide() == 2) {OtherSide = 1;}

	if (fArray[cal->GetDetector()-1][cal->GetFace()-1][OtherSide-1][cal->GetStrip()-1] != -1) {index2 = 1;}
	for (Int_t jj = 0; jj < fNumStrips; jj++)
	{
		if (fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][jj] != -1 && fArray[cal->GetDetector()-1][OtherFace-1][OtherSide-1][jj] != -1) 
		{
			index3 = 1; index4 = 1; break;
		}
	}
	//printf("I exist %d %d %d\n",index2,index3,index4);
	

        //printf("Det: %d, Face: %d, Side: %d, Strip: %d, Index: %d\n",cal->GetDetector()-1,OtherFace-1,OtherSide-1,cal->GetStrip(), fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][cal->GetStrip()-1]);
	if (index2 == -1 || index3 == -1 || index4 == -1) {return kFALSE;}
	else return kTRUE;
}

Int_t R3BPspxCal2Hit::FacePartner(Int_t index) 
{
	auto cal = (R3BPspxCalData const*)fCal->At(index);

	Int_t Time=-1000; Int_t Time_partner=-9000;

	Int_t OtherFace = -1;
       	Int_t OtherSide = -1;
	if (cal->GetFace() == 1) {OtherFace = 2;}
	else if (cal->GetFace() == 2) {OtherFace = 1;}
	if (cal->GetSide() == 1) {OtherSide = 2;}
	else if (cal->GetSide() == 2) {OtherSide = 1;}
	//static int faceindex[MultFace(cal->GetDetector()-1,OtherFace-1)];
	Int_t faceindex=-1;


	for (Int_t jj = 0; jj < fNumStrips; jj++)
	{
		if (fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][jj] != -1 && fArray[cal->GetDetector()-1][OtherFace-1][OtherSide-1][jj] != -1) 
		{
			auto calpartner = (R3BPspxCalData const*)fCal->At(fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][jj]);
			auto calpartner2 = (R3BPspxCalData const*)fCal->At(fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][jj]);

			Int_t TotalEnergy = calpartner->GetEnergy() + calpartner2->GetEnergy();
			if (TotalEnergy > 0.1*maxenergy[cal->GetDetector()-1]) {

				Time_partner = GetTimeCorr_s473(calpartner->GetTime()); 
				Time = GetTimeCorr_s473(cal->GetTime());
				if ((Time - Time_partner) < 15 && (Time - Time_partner) > -15)
				{	
				
					faceindex = fArray[cal->GetDetector()-1][OtherFace-1][cal->GetSide()-1][jj]; 
				}
			}
		}
	}
	return faceindex;
}

Int_t R3BPspxCal2Hit::MultFace(Int_t det, Int_t face)
{
	Int_t mult = 0;

	for (Int_t ii = 0; ii < fNumStrips; ii++)
	{
		if (fArray[det][face][1][ii] != -1) 
		{
			if (MatchedEvents(fArray[det][face][1][ii]) == kTRUE) {mult++;}
		}
	}
	return mult;
}

void R3BPspxCal2Hit::SetNumExpt(Int_t expt) {fNumExpt = expt;}
void R3BPspxCal2Hit::SetNumDet(Int_t det) {fNumDet = det;}
void R3BPspxCal2Hit::SetNumStrips(Int_t strip) {fNumStrips = strip;}

ClassImp(R3BPspxCal2Hit)
