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

// ------------------------------------------------------------------------
// -----                    R3BPspxCal2HitPar                      -----
// -----              Created  Jul 2022 by A. Jedele                  -----
// -----      Adopted from script by M. Holl, S. Storck                 -----
// -----  Purpose of script is to generate the position calibration   -----
// -----        and the timing calibration (new addition)             -----
// ------------------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "FairRunAna.h"
#include "FairRootManager.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxCal2HitPar.h"
#include "R3BPspxCalData.h"
#include "R3BPspxHitPar.h"

#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <stdio.h>

//ROOT headers
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TDirectory.h"
#include "TProfile.h"

//R3BPspsxCal2HitPar: Default Constructor ---------------------------------------
R3BPspxCal2HitPar::R3BPspxCal2HitPar()
    : FairTask("PspxCal2HitPar", 1)
{
}

//R3BPspsxCal2HitPar: Standard Constructor ---------------------------------------
R3BPspxCal2HitPar::R3BPspxCal2HitPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCal(nullptr)
    , fParOutName("DefaultOutput")
    , fNumDet(0)
    , fNumStrips(0)
{
}

//Virtual R3BPspsxCal2HitPar: Destructor ----------------------------------------
R3BPspxCal2HitPar::~R3BPspxCal2HitPar()
{
    LOG(INFO) << "R3BPspxCal2HitPar: Delete instance";
    delete fCal;
}

//----Public method Init -----------------------------------------------------------
InitStatus R3BPspxCal2HitPar::Init()
{
	LOG(INFO) << "R3BPspxCal2HitPar: Init" ;
	
	//FairRootManager
	FairRootManager* fMan = FairRootManager::Instance();
	if (!fMan) { LOG(ERROR) << "R3BPspxCal2HitPar::Init(). No FairRootManager found"; return kFATAL;}
	
	//Get Cal Data
	 fCal = (TClonesArray*)fMan->GetObject("PspxCalData");
	if (!fCal) {LOG(ERROR) << "R3BPspxCal2HitPar::Init(). No PspxCalData found."; return kFATAL;}
	
	//Get Cal Data
	
	// R3BEventHeader for trigger information, if needed!
	//fHeader = (R3BEventHeader*)fMan->GetObject("R3BEventHeader");
	
	//Container needs to be created in tcal/R3BCalContFact.cxx AND R3BCal needs
	//to be set as dependency in CMakelists.txt (in this case in the psp dir.)
	FairRuntimeDb* rtdb = FairRuntimeDb::instance();
	if (!rtdb) {return kFATAL;}
	
	if (!fNumDet)
	{
	    LOG(ERROR) << "R3BPspxCal2HitPar::Init(). No Detectors detected.";
	    return kFATAL;
	}
	
	if (!fNumStrips)
	{
	    LOG(ERROR) << "R3BPspxCal2HitPar::Init(). No Strips found.";
	    return kFATAL;
	}

	outputfile.open(fParOutName.Data());
	
	//Define TGraph for the fits and other info
	//4 Parameters per strip per detectors: 2 sides and 2 faces
	printf("\n\nNumDet: %d, NumStrips: %d\n\n\n",fNumDet,fNumStrips);
	TString Name;
	for (Int_t ii = 0; ii < fNumDet; ii++)
	{
		for  (Int_t jj = 0; jj<fNumStrips; jj++)
		{
			for  (Int_t kk = 0; kk<fNumStrips; kk++)
			{
				energydir->
				Name = Form("RawEnergyFront_Det%d_Strip1%d_Strip2%d",ii+1,jj+1,kk+1);
				RawEnergyFront[ii][jj][kk] = new TH1F(Name.Data(),Name.Data(),4000,0,4000000);
				Name = Form("RawEnergyBack_Det%d_Strip1%d_Strip2%d",ii+1,jj+1,kk+1);
				RawEnergyBack[ii][jj][kk] = new TH1F(Name.Data(),Name.Data(),4000,0,4000000);
			
				//if ((ii%2==0) && ())
				//{
					//Name = Form("RawEnergyTime_Front_Det%d_Strip1%dStrip2%d",ii+1,jj+1,kk+1);
					//RawEnergyTimeFront[ii][jj][kk] = new TH2D(Name.Data(),Name.Data(),2048,0,4095,4000,0,4000000);
					//Name = Form("RawEnergyTime_Back_Det%d_Strip1%dStrip2%d",ii+1,jj+1,kk+1);
					//RawEnergyTimeBack[ii][jj][kk] = new TH2D(Name.Data(),Name.Data(),2048,0,4095,4000,0,4000000);
				//}
			}
			Name = Form("RawEnergyFront_Det%d_Strip%d",ii+1,jj+1);
			RawEnergyFront2D[ii][jj] = new TH2F(Name.Data(),Name.Data(),32,1,32,4000,0,4000000);
			Name = Form("RawEnergyBack_Det%d_Strip%d",ii+1,jj+1);
			RawEnergyBack2D[ii][jj] = new TH2F(Name.Data(),Name.Data(),32,1,32,4000,0,4000000);
		}
	}
	return kSUCCESS;
}

//----Public method ReInit -----------------------------------------------------------
InitStatus R3BPspxCal2HitPar::ReInit() { return kSUCCESS; }

void R3BPspxCal2HitPar::Exec(Option_t* option)
{
	Reset();
    /**
     * Gets the parameters for the conversion from Cal to Hit level. It is called for every event.
     * The parameters are energy and timing.
     * Energy: 2 step calibration: 1st step is to align all the corrected energies. 2nd is to find the offset and slope using various different energies
     * Timing: calibration was done in the 1st step
     * Sets beam peak to 3x10^6
     *
     * Face 1 is strips in the x-direction (resitive splitting gives y-position)
     * Face 2 is strips in the y-direction (resitive splitting gives x-position)
     *
     * Side 1 is left, Side 2 is right (look beam-downstream)
     * Side 1 is up, Side 2 is down (look beam-downstream)
     *
     * Face 1 is biased
     */

	counter_tot++;

	Int_t nHits = fCal->GetEntries();
	if (nHits > 0) 
	{
	       	counter_events++; /*PRINTF("entries: %d\n",fCal->GetEntries());*/

		//Start of the actual calibration
		for (Int_t ii=0; ii<nHits; ii++)
		{
        		auto cal = (R3BPspxCalData const*)fCal->At(ii);

			Float_t TotalEnergy_Front = cal->GetEnergyFront1() + cal->GetEnergyFront2();
			Float_t TotalEnergy_Back  = cal->GetEnergyBack1() + cal->GetEnergyBack2();
			RawEnergyFront[cal->GetDetector()-1][cal->GetStripFront()-1][cal->GetStripBack()-1]->Fill(TotalEnergy_Front);
			RawEnergyBack[cal->GetDetector()-1][cal->GetStripBack()-1][cal->GetStripFront()-1]->Fill(TotalEnergy_Back);
			RawEnergyFront2D[cal->GetDetector()-1][cal->GetStripFront()-1]->Fill(cal->GetStripBack()-1,TotalEnergy_Front);
			RawEnergyBack2D[cal->GetDetector()-1][cal->GetStripBack()-1]->Fill(cal->GetStripFront()-1,TotalEnergy_Back);
			//printf("Testing nHits. Det: %d, StripFront: %d, StripBack: %d, TotalEnergy_Front: %lf, TotalEnergy_Back: %lf\n",cal->GetDetector(),cal->GetStripFront(),cal->GetStripBack(),TotalEnergy_Front, TotalEnergy_Back);
		}
	}
	fCal->Clear();
}


void R3BPspxCal2HitPar::Reset() {}


void R3BPspxCal2HitPar::FinishTask() 
{
	printf("R3BPspxCal2Hit::FinishTask()\n");	
	GetEnergyCorrection();
	for (Int_t ii=0; ii < fNumDet; ii++)
	{
		for (Int_t jj=0; jj < fNumStrips; jj++)
		{
			for (Int_t ll=0; ll < fNumStrips; ll++)
			{
				if(RawEnergyFront[ii][jj][ll]->GetEntries()>0) {RawEnergyFront[ii][jj][ll]->Write();}	
				if(RawEnergyBack[ii][jj][ll]->GetEntries()>0)  {RawEnergyBack[ii][jj][ll]->Write();}
				fArray[ii][jj][ll]=-1;	
			}
			if(RawEnergyFront2D[ii][jj]->GetEntries() > 0) RawEnergyFront2D[ii][jj]->Write();
			if(RawEnergyBack2D[ii][jj]->GetEntries() > 0) RawEnergyBack2D[ii][jj]->Write();
		}
	}
	outputfile.close();
	printf("Counters. Total: %d, Events: %d\n", counter_tot, counter_events);
}

void R3BPspxCal2HitPar::GetEnergyCorrection() 
{
	printf("Inside GetEnergyCorrection\n");
	for (int ii = 0; ii < fNumDet; ii++)
	{
		TString nom = Form("GausFits%d",ii);
		GausFitCan[ii] = new TCanvas(nom.Data(),nom.Data(),800,800);
		Int_t graphcounter = 0;
		
		for (int jj = 0; jj < fNumStrips; jj++)
		{
			TGraphErrors *graph = new TGraphErrors();
			Int_t counter = 0;
			for (int kk = 0; kk < fNumStrips; kk++)
			{
				//Float_t centermax = 0.; Float_t CenterPar = 0.;
				//if (RawEnergyFront[ii][15][kk]->GetEntries() > 20)
				//{
				//	centermax = RawEnergyFront[ii][15][kk]->GetMaximumBin();
				//	TF1 *GausFitCenter = new TF1("Gaussian Fit Center","gaus(0)",0,6000000);
				//	RawEnergyFront[ii][15][kk]->Fit(GausFitCenter,"QNEW0","",(centermax*10-1000),(centermax*10+1000));
				//	CenterPar = (Float_t)GausFitCenter->GetParameter(1);
				//}
				if (RawEnergyFront[ii][jj][kk]->GetEntries() > 50)
				{
					Float_t max = RawEnergyFront[ii][jj][kk]->GetMaximumBin();
					printf("Max bins: %d\n", RawEnergyFront[ii][jj][kk]->GetMaximumBin());
					TF1 *GausFit = new TF1("Gaussian Fit","gaus(0)",0,4000000);
					if (ii%2 == 0) {RawEnergyFront[ii][jj][kk]->Fit(GausFit,"QNEW0","",1000*(max-10),1000*(max+10));}
					else RawEnergyFront[ii][jj][kk]->Fit(GausFit,"QNEW0","",1000*(max-10),1000*(max+10));

					graph->SetPoint(counter,kk+1,(Float_t)GausFit->GetParameter(1));

					//fArray[ii][jj][kk] = CenterPar/(Float_t)GausFit->GetParameter(1);
					printf("det: %d, StripFront: %d, StripBack: %d, Gaussian mean: %lf\n",ii+1,jj+1,kk+1,GausFit->GetParameter(1));
					//outputfile << ii << "\t" << 1 << "\t" << jj << "\t" << (Int_t)(XX+0.25) << "\t" << ratioY << std::endl;
					counter++;
				}
				//centermax = 0.;
			}
			//TF1 *QuadFit = new TF1("Quadratic Fit","[0]+[1]*(x)+[2]*(x)*(x)",0,32);
			//QuadFit->SetParLimits();
			//if (graph->GetN()>2) graph->Fit(QuadFit,"NEW0");
			//QuadFit->Draw();
			if (graph->GetN() > 0) 
			{
				if (graphcounter == 0) {graph->Draw("AP");}
				else {graph->Draw("P");}
				graph->SetMarkerColor(jj+1);
				graph->SetMarkerStyle(20);
				graphcounter++;
			}
		}
	}
};

void R3BPspxCal2HitPar::SetNumDet(Int_t det) {fNumDet = det;}
void R3BPspxCal2HitPar::SetNumStrips(Int_t strip) {fNumStrips = strip;}
void R3BPspxCal2HitPar::SetParOutName(TString paroutname) {fParOutName = paroutname;}

ClassImp(R3BPspxCal2HitPar)
