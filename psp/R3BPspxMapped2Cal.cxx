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
// -----                 R3BPspxMapped2Cal              -----
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
#include "R3BPspxMapped2Cal.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxCalData.h"
#include "R3BPspxCalPar.h"

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

R3BPspxMapped2Cal::R3BPspxMapped2Cal()
    : FairTask("PspxMapped2Cal", 1)
{
}

R3BPspxMapped2Cal::R3BPspxMapped2Cal(const char* name, Int_t iVerbose)
    	: FairTask(name, iVerbose)
    	, fNumExpt(0)
    	, fNumDet(0)
	, fNumStrips(0)
//	, fCalPar(nullptr)
	, fMapped(nullptr)
	, fCal(nullptr)
{
}

R3BPspxMapped2Cal::~R3BPspxMapped2Cal()
{
//	delete fCalPar;
	delete fMapped;
	delete fCal;
}

void R3BPspxMapped2Cal::SetParContainers()
{
	//Parameter Container
	//Reading PspCalPar from FairRuntimeDB
	FairRuntimeDb* rtdb = FairRuntimeDb::instance();
	if (!rtdb)
	{
		LOG(ERROR) << "FairRuntimeDb not found";
	}

	//fCalPar = (R3BPspxCalPar*)rtdb->getContainer("PspxCalPar");
	//if (!fCalPar)
	//{
	//	LOG(ERROR) << "R3BPspxMapped2Cal::SetParContainers() Couldn't get handle on pspxCalPar container";
	//}
	//else	
	//{
	//	LOG(INFO) << "R3BPspxMapped2Cal::SetParContainers() pspxCalPar container open";
	//}
}

InitStatus R3BPspxMapped2Cal::Init()
{
	//FairRootManager
	FairRootManager* fMan = FairRootManager::Instance();
	if (!fMan) { LOG(ERROR) << "R3BPspxMapped2Cal::Init(). No FairRootManager found"; return kFATAL;}
	
	//Get Mapped Data
	fMapped = (TClonesArray*)fMan->GetObject("PspxMapped");
	if (!fMapped) {LOG(ERROR) << "R3BPspxMapped2Cal::Init(). No PspxMappedData found."; return kFATAL;}
	printf("fMapped->GetEntries(): %d\n", fMapped->GetEntries());

	// R3BEventHeader for trigger information, if needed!
	//fHeader = (R3BEventHeader*)fMan->GetObject("R3BEventHeader");
	
	//Container needs to be created in tcal/R3BCalContFact.cxx AND R3BCal needs
	//to be set as dependency in CMakelists.txt (in this case in the psp dir.)
	FairRuntimeDb* rtdb = FairRuntimeDb::instance();
	if (!rtdb) {return kFATAL;}
	
//	//Get the Cal Par container
//	fCalPar = (R3BPspxCalPar*)rtdb->getContainer("PspxCalPar");
//	if (!fCalPar) {
//	    LOG(ERROR) << "R3BPspxMapped2Cal::Init(). Could get handle of PspxCalPar ";
//	    return kFATAL;
//	}
	
	if (!fNumDet)
	{
	    LOG(ERROR) << "R3BPspxMapped2Cal::Init(). No Detectors detected.";
	    return kFATAL;
	}
	
	if (!fNumStrips)
	{
	    LOG(ERROR) << "R3BPspxMapped2Cal::Init(). No Strips found.";
	    return kFATAL;
	}


	//Register Cal Level Data
	fCal = new TClonesArray("R3BPspxCalData");
	fMan->Register("PspxCalData","PSPX Cal", fCal, kTRUE);
	fCal->Clear();
	
	printf("\n\nNumDet: %d, NumStrips: %d\n\n\n",fNumDet,fNumStrips);
	TString Name;
	for (Int_t ii = 0; ii < fNumDet; ii++)
	{
    		for (Int_t jj = 0; jj < 2; jj++)
    		{
	   		for  (Int_t ll = 0; ll<fNumStrips; ll++)
    			{
           			Name = Form("TimeFront_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hTimeFront[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(),32,1,33,4096,0,4096);
				Name = Form("TimeBack_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hTimeBack[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(),32,1,33,4096,0,4096);
 
           			Name = Form("EnergyFront_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hEnergyFront[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(),32,1,33,10000,0,1000000);
 
           			Name = Form("EnergyCorrFront_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hEnergyCorrFront[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(), 32,1,33,10000,0,1000000);
           			
				Name = Form("EnergyBack_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hEnergyBack[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(),32,1,33,10000,0,1000000);
 
           			Name = Form("EnergyCorrBack_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    			hEnergyCorrBack[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(), 32,1,33,10000,0,1000000);
				for (Int_t kk = 0; kk < fNumStrips; kk++)
				{
					GainFactors[ii][jj][kk][ll] = -1.;
				}
			}
		}
	}

	TString filename = Form("Mapped2CalParGlobal.txt");
	std::ifstream infile(filename);
	if (infile.bad()) {printf("Zombie File!!!\n");}

	Int_t Det=0, Face=0, StripFront=0, StripBack=0;
	Double_t ratio=0.;
	Char_t line[1000];

	while(!infile.eof())
	{
		infile.getline(line,1000);
		sscanf(line,"%d %d %d %d %lf", &Det, &Face, &StripFront, &StripBack, &ratio);
		GainFactors[Det][Face][StripFront][StripBack] = ratio;
		//if (Det != 0) {printf("SScanf line. Det; %d, Face: %d, StripFront: %d, StripBack: %d, ratio: %lf, GainFactor: %lf\n", Det, Face, StripFront, StripBack, ratio, GainFactors[Det][Face][StripFront][StripBack]);}
	}


    return kSUCCESS;
}

InitStatus R3BPspxMapped2Cal::ReInit() {return kSUCCESS;}

void R3BPspxMapped2Cal::Exec(Option_t* option)
{
	Int_t nHits = fMapped->GetEntries();
	if (nHits > 0)
	{
		for (Int_t ii = 0; ii < nHits; ii++)
    		{
			auto mapped = (R3BPspxMappedData const*)fMapped->At(ii);
			if (mapped == NULL) {ZombieAreAlive = kTRUE; printf("I'm am Zombie\n"); break;}  // removes bad events where -nan is filled into the fMapped array
			fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1] = ii; 
		}

		for (Int_t ii=0; ii<nHits; ii++)
		{
			if (ZombieAreAlive == kTRUE) {break;}

			counter_total++;
			auto mapped = (R3BPspxMappedData const*)fMapped->At(ii);
			if (MatchedEvents(ii) == kFALSE) {continue;}

			//Used to get the array index for the opposite side and opposite face of a detector
			Int_t OtherSide = 0; Int_t OtherFace = 0;
			if (mapped->GetSide() == 1) {OtherSide = 2;}
			else if (mapped->GetSide() == 2) {OtherSide = 1;}
			if (mapped->GetFace() == 1) {OtherFace = 2;}
			else if (mapped->GetFace() == 2) {OtherFace = 1;}
			auto mapped2 = (R3BPspxMappedData const*)fMapped->At(fArray[mapped->GetDetector()-1][mapped->GetFace()-1][OtherSide-1][mapped->GetStrip()-1]);	

			Float_t time_corr1=0; Float_t time_corr2=0;
			Float_t time_corr3=0; Float_t time_corr4=0;

			//Code specific for Nik's FW
			if (fNumExpt == 473 || fNumExpt == 444)
			{
				time_corr1 = GetTimeCorr_s473(mapped->GetTime());
				time_corr2 = GetTimeCorr_s473(mapped2->GetTime());
			}

			//Code specific for CALIFA FW 
			if (fNumExpt == 515)
			{
				time_corr1 = GetTimeCorr_s515(mapped->GetTime(),mapped->GetTrigger());
				time_corr2 = GetTimeCorr_s515(mapped2->GetTime(),mapped->GetTrigger());
			}

			Float_t timediff = time_corr1 - time_corr2;

			if (timediff > -10 && timediff < 10)
			{

				//Multiplicity requirement: Mult x-face has to equal mult y-face - timing and energy calibration needs this requirement
				if ((MultFace(mapped->GetDetector()-1,mapped->GetFace()-1) == 1) && (MultFace(mapped->GetDetector()-1,OtherFace-1) == 1)) 
				{
					if (mapped->GetSide()==1 && mapped->GetFace()==1)
					{
						if (FacePartner(ii) == -1) {continue;}
						else 
						{
							auto mappedback = (R3BPspxMappedData const*)fMapped->At(FacePartner(ii));
							Int_t OtherSideBack = 0;
							if (mappedback->GetSide() == 1) {OtherSideBack = 2;}
							else if (mappedback->GetSide() == 2) {OtherSideBack = 1;}
							auto mappedback2 = (R3BPspxMappedData const*)fMapped->At(fArray[mappedback->GetDetector()-1][mappedback->GetFace()-1][OtherSideBack-1][mappedback->GetStrip()-1]);
							//Code specific for Nik's FW
							if (fNumExpt == 473 || fNumExpt == 444)
							{
								time_corr3 = GetTimeCorr_s473(mappedback->GetTime());
								time_corr4 = GetTimeCorr_s473(mappedback2->GetTime());
							}

							//Code specific for CALIFA FW 
							if (fNumExpt == 515)
							{
								time_corr3 = GetTimeCorr_s515(mappedback->GetTime(),mapped->GetTrigger());
								time_corr4 = GetTimeCorr_s515(mappedback2->GetTime(),mapped->GetTrigger());
							}
							if (GainFactors[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1][mappedback->GetStrip()-1] != -1 && GainFactors[mappedback->GetDetector()-1][mappedback->GetFace()-1][mappedback->GetStrip()-1][mapped->GetStrip()-1] != -1)
							{
								//printf("Det: %d, face: %d, FrontStrip: %d, BackStrip: %d\n",mapped->GetDetector()-1,mapped->GetFace()-1,mapped->GetStrip()-1,mappedback->GetStrip()-1);
								hTimeFront[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mappedback->GetStrip()-1, mapped->GetTime());
								hTimeBack[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]-> Fill(mapped->GetStrip()-1, mappedback->GetTime());
								hEnergyFront[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mappedback->GetStrip()-1, mapped->GetEnergy() + mapped2->GetEnergy());
								hEnergyCorrFront[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mappedback->GetStrip()-1, mapped->GetEnergy() + mapped2->GetEnergy() * GainFactors[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1][mappedback->GetStrip()-1]);
								hEnergyBack[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mappedback->GetStrip()-1, mapped->GetEnergy() + mapped2->GetEnergy());
								hEnergyCorrBack[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mappedback->GetStrip()-1, mapped->GetEnergy() + mapped2->GetEnergy() * GainFactors[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1][mappedback->GetStrip()-1]);
								//printf("Det: %d, Face: %d, Strip: %d, Back Strip: %d, Energy1: %d, Energy2: %d, GainFactor: %lf\n", mapped->GetDetector()-1, mapped->GetFace()-1, mapped->GetStrip()-1, mappedback->GetStrip()-1, mapped->GetEnergy(), mapped2->GetEnergy(), GainFactors[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1][mappedback->GetStrip()-1]);
        							//printf("Det: %d, Face: %d, Side: %d, Strip: %d, Index: %d, nHits: %d\n",mapped->GetDetector()-1,OtherFace-1,OtherSide-1,mapped->GetStrip(), fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1],nHits);
								new ((*fCal)[fCal->GetEntriesFast()]) R3BPspxCalData(mapped->GetDetector(), mapped->GetStrip(), mappedback->GetStrip(), (Float_t)mapped->GetEnergy(), (Float_t)mapped2->GetEnergy()*(Float_t)GainFactors[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1][mappedback->GetStrip()-1], (Float_t)mappedback->GetEnergy(), (Float_t)mappedback2->GetEnergy()*(Float_t)GainFactors[mappedback->GetDetector()-1][mappedback->GetFace()-1][mappedback->GetStrip()-1][mapped->GetStrip()-1], time_corr1, time_corr2, time_corr3, time_corr4, MultFace(mapped->GetDetector()-1,mapped->GetFace()-1), MultFace(mappedback->GetDetector()-1,mappedback->GetFace()-1));
							
								//if (fCal->GetEntriesFast() % 1000 == 0) printf("Entries: %d, Counter: %d, Counter_Total: %d, weird one: %d\n",fCal->GetEntriesFast(), counter, counter_total, counter*fCal->GetEntries());
								counter++;
							}
							else continue;
						}
					}
				}
			}
		}
	}
}

void R3BPspxMapped2Cal::FinishEvent() 
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
	fCal->Clear();
}

void R3BPspxMapped2Cal::FinishTask()
{
	for (Int_t ii=0; ii < fNumDet; ii++)
	{
		for (Int_t jj=0; jj < 2; jj++)
		{
			for (Int_t ll=0; ll < fNumStrips; ll++)
			{
				hEnergyFront[ii][jj][ll]->Write();
				hEnergyCorrFront[ii][jj][ll]->Write();
				hEnergyBack[ii][jj][ll]->Write();
				hEnergyCorrBack[ii][jj][ll]->Write();
			}
		}
	}
}

Float_t R3BPspxMapped2Cal::GetTimeCorr_s473(Int_t time)
{
	//corrects for the timing offset
	Int_t corr_time = time - 2048;
	if (corr_time < 0) {corr_time = 4096 + corr_time;}
	return corr_time;
}

Float_t R3BPspxMapped2Cal::GetTimeCorr_s515(Int_t time, Int_t trigger)
{
	//corrects for the timing offset
	Int_t corr_time = time - trigger;
	return corr_time;
}

Bool_t R3BPspxMapped2Cal::MatchedEvents(Int_t index) 
{
	Int_t index2 = -1, index3 = -1, index4 = -1;
	auto mapped = (R3BPspxMappedData const*)fMapped->At(index);

        //printf("Event Number: %d, Det: %d, Face: %d, Side: %d, Strip: %d, Index: %d\n",counter_tot,mapped->GetDetector(),mapped->GetFace(),mapped->GetSide(),mapped->GetStrip(), fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1]);
	if (fMapped->GetEntries() < 4) {return kFALSE;}
	Int_t OtherFace = -1;
	Int_t OtherSide = -1;
	if (mapped->GetFace() == 1) {OtherFace = 2;}
	else if (mapped->GetFace() == 2) {OtherFace = 1;}
	if (mapped->GetSide() == 1) {OtherSide = 2;}
	else if (mapped->GetSide() == 2) {OtherSide = 1;}

	if (fArray[mapped->GetDetector()-1][mapped->GetFace()-1][OtherSide-1][mapped->GetStrip()-1] != -1) {index2 = 1;}
	for (Int_t jj = 0; jj < fNumStrips; jj++)
	{
		if (fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj] != -1 && fArray[mapped->GetDetector()-1][OtherFace-1][OtherSide-1][jj] != -1) 
		{
			index3 = 1; index4 = 1; break;
		}
	}
	//printf("I exist %d %d %d\n",index2,index3,index4);
	

        //printf("Det: %d, Face: %d, Side: %d, Strip: %d, Index: %d\n",mapped->GetDetector()-1,OtherFace-1,OtherSide-1,mapped->GetStrip(), fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][mapped->GetStrip()-1]);
	if (index2 == -1 || index3 == -1 || index4 == -1) {return kFALSE;}
	else return kTRUE;
}

Int_t R3BPspxMapped2Cal::FacePartner(Int_t index) 
{
	auto mapped = (R3BPspxMappedData const*)fMapped->At(index);

	Int_t Time=-1000; Int_t Time_partner=-9000;

	Int_t OtherFace = -1;
       	Int_t OtherSide = -1;
	if (mapped->GetFace() == 1) {OtherFace = 2;}
	else if (mapped->GetFace() == 2) {OtherFace = 1;}
	if (mapped->GetSide() == 1) {OtherSide = 2;}
	else if (mapped->GetSide() == 2) {OtherSide = 1;}
	//static int faceindex[MultFace(mapped->GetDetector()-1,OtherFace-1)];
	Int_t faceindex=-1;


	for (Int_t jj = 0; jj < fNumStrips; jj++)
	{
		if (fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj] != -1 && fArray[mapped->GetDetector()-1][OtherFace-1][OtherSide-1][jj] != -1) 
		{
			auto mappedpartner = (R3BPspxMappedData const*)fMapped->At(fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj]);
			auto mappedpartner2 = (R3BPspxMappedData const*)fMapped->At(fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj]);

			Int_t TotalEnergy = mappedpartner->GetEnergy() + mappedpartner2->GetEnergy();
			if (TotalEnergy > 0.1*maxenergy[mapped->GetDetector()-1]) {

				Time_partner = GetTimeCorr_s473(mappedpartner->GetTime()); 
				Time = GetTimeCorr_s473(mapped->GetTime());
				if ((Time - Time_partner) < 15 && (Time - Time_partner) > -15)
				{	
				
					faceindex = fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj]; 
				}
			}
		}
	}
	return faceindex;
}

Int_t R3BPspxMapped2Cal::MultFace(Int_t det, Int_t face)
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

void R3BPspxMapped2Cal::SetNumExpt(Int_t expt) {fNumExpt = expt;}
void R3BPspxMapped2Cal::SetNumDet(Int_t det) {fNumDet = det;}
void R3BPspxMapped2Cal::SetNumStrips(Int_t strip) {fNumStrips = strip;}

ClassImp(R3BPspxMapped2Cal)
