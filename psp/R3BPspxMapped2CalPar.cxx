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
// -----                    R3BPspxMapped2CalPar                      -----
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
#include "R3BPspxMapped2CalPar.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxCalPar.h"
#include "R3BPspxContFact.h"

#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <stdio.h>

//ROOT headers
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TDirectory.h"
#include "TProfile.h"

//R3BPspsxMapped2CalPar: Default Constructor ---------------------------------------
R3BPspxMapped2CalPar::R3BPspxMapped2CalPar()
    : FairTask("PspxMapped2CalPar", 1)
{
}

//R3BPspsxMapped2CalPar: Standard Constructor ---------------------------------------
R3BPspxMapped2CalPar::R3BPspxMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapped(nullptr)
    , fParOutName("DefaultOutput")
    , fNumDet(0)
    , fNumStrips(0)
    , fNumExpt(0)
{
}

//Virtual R3BPspsxMapped2CalPar: Destructor ----------------------------------------
R3BPspxMapped2CalPar::~R3BPspxMapped2CalPar()
{
    LOG(INFO) << "R3BPspxMapped2CalPar: Delete instance";
    delete fMapped;
}

//----Public method Init -----------------------------------------------------------
InitStatus R3BPspxMapped2CalPar::Init()
{
	LOG(INFO) << "R3BPspxMapped2CalPar: Init" ;
	
	//FairRootManager
	FairRootManager* fMan = FairRootManager::Instance();
	if (!fMan) { LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No FairRootManager found"; return kFATAL;}
	
	//Get Mapped Data
	fMapped = (TClonesArray*)fMan->GetObject("PspxMapped");
	if (!fMapped) {LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No PspxMappedData found."; return kFATAL;}
	
	//Container needs to be created in tcal/R3BCalContFact.cxx AND R3BCal needs
	//to be set as dependency in CMakelists.txt (in this case in the psp dir.)
	FairRuntimeDb* rtdb = FairRuntimeDb::instance();
	if (!rtdb) {return kFATAL;}
	
	if (!fNumDet)
	{
	    LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No Detectors detected.";
	    return kFATAL;
	}
	
	if (!fNumStrips)
	{
	    LOG(ERROR) << "R3BPspxMapped2CalPar::Init(). No Strips found.";
	    return kFATAL;
	}
	

	outputfile.open(fParOutName.Data());

	//Define TGraph for the fits and other info
	//4 Parameters per strip per detectors: 2 sides and 2 faces
	printf("\n\nNumDet: %d, NumStrips: %d\n\n\n",fNumDet,fNumStrips);
	TString Name;
	for (Int_t ii = 0; ii < fNumDet; ii++)
	{
    		for (Int_t jj = 0; jj < 2; jj++)
    		{
    			for (Int_t kk = 0; kk < 2; kk++)
    			{
	   			for  (Int_t ll = 0; ll<fNumStrips; ll++)
    				{
           				Name = Form("Time_Det%d_Face%d_Side%d_Strip%d",ii+1,jj+1,kk+1,ll+1);
	    				hTime[ii][jj][kk][ll] = new TH1F(Name.Data(), Name.Data(), 5000,-5000,5000);
					Energy->cd();
           				Name = Form("Energy_Det%d_Face%d_Side%d_Strip%d",ii+1,jj+1,kk+1,ll+1);
	    				hEnergy[ii][jj][kk][ll] = new TH1F(Name.Data(), Name.Data(), 10000,0,1000000);
 
           				if (kk == 0)
					{
						Energy->cd();
           					Name = Form("EnergyTot_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    					hEnergyTot[ii][jj][ll] = new TH1F(Name.Data(), Name.Data(), 10000,0,1000000);
 
						Interstrip->cd();
           					Name = Form("Interstrip_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    					hInterstrip[ii][jj][ll] = new TH1F(Name.Data(), Name.Data(), 10000,0,1000000);
           					Name = Form("InterstripTot_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    					hInterstripTot[ii][jj][ll] = new TH1F(Name.Data(), Name.Data(), 10000,0,1000000);
           					Name = Form("Strip1vsStrip2_Det%d_Face%d_Strip%d",ii+1,jj+1,ll+1);
	    					hStrip1vsStrip2[ii][jj][ll] = new TH2F(Name.Data(), Name.Data(),1000,-1.0,1.0,1000,-1.0,1.0);
	   
	
						Position->cd();
       						Name = Form("PosX_Det%d_Strip%d",ii+1,2*ll+jj+1);
						hPosX_strip[ii][2*ll+jj]   = new TH1F(Name.Data(), Name.Data(),1000,-1.0,1.0);
        					Name = Form("PosY_Det%d_Strip%d",ii+1,2*ll+jj+1);
						hPosY_strip[ii][2*ll+jj]   = new TH1F(Name.Data(), Name.Data(),1000,-1.0,1.0);

						if (jj==0)
						{
							Name = Form("PosX_Det%dStrip%d",ii+1,ll+1);
							hPosX[ii][ll] = new TH2F(Name.Data(),Name.Data(),64,1.0,33.0,1000,-1.0,1.0);
							Name = Form("PosY_Det%dStrip%d",ii+1,ll+1);
							hPosY[ii][ll] = new TH2F(Name.Data(),Name.Data(),64,1.,33.,1000,-1.0,1.0);
						}
					}
				}
			}
	    	}
		Position->cd();
		Name = Form("PosXvsStrip_Det%d",ii+1);
		hPosXvsStrip_det[ii]   = new TH2F(Name.Data(), Name.Data(),64,1.0,33.0,1000,-1.0,1.0);
		Name = Form("PosYvsStrip_Det%d",ii+1);
		hPosYvsStrip_det[ii]   = new TH2F(Name.Data(), Name.Data(),64,1.0,33.0,1000,-1.0,1.0);
		
		Name = Form("PosX_Det%d",ii+1);
		hPosX_det[ii] = new TH2F(Name.Data(),Name.Data(),1000,-1.0,1.0,64,1.0,33.0);
		Name = Form("PosY_Det%d",ii+1);
		hPosY_det[ii] = new TH2F(Name.Data(),Name.Data(),64,1.,33.,1000,-1.0,1.0);
		
		Name = Form("HitMap_Det%d",ii+1);
		hHitMap[ii] = new TH2F(Name.Data(),Name.Data(),32,1,33,32,1,33);
		
		Name = Form("TimeDiffPartner_Det%d",ii+1);
		hTimeDiffPartner[ii] = new TH1F(Name.Data(),Name.Data(),1000,-500,500);
		
		Name = Form("Mult_Det%d",ii+1);
		hMult[ii] = new TH2F(Name.Data(),Name.Data(),100,0,100,100,0,100);
	}	
	return kSUCCESS;
}

//----Public method ReInit -----------------------------------------------------------
InitStatus R3BPspxMapped2CalPar::ReInit() { return kSUCCESS; }

void R3BPspxMapped2CalPar::Exec(Option_t* option)
{
    /**
     * Gets the parameters for the conversion from Mapped to Cal level. It is called for every event.
     * The parameters are energy, timing and position.
     * Energy: thresholds found and applied
     * Timing: determines if the timing difference between each side of the strip is approx. the same
     * Sets beam peak to 0
     * Position: use interstrip events to determine the position stretching parameters
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

	Int_t nHits = fMapped->GetEntries();
	if (nHits > 0) 
	{
	       	counter_events++; /*PRINTF("entries: %d\n",fMapped->GetEntries());*/

		//Fill index Array for all hits - used to calculate multiplicity, find partners, identify interstrip events
		for (Int_t ii = 0; ii < nHits; ii++)
		{
			auto mapped = (R3BPspxMappedData const*)fMapped->At(ii);
			if (mapped == NULL) {ZombieAreAlive = kTRUE; printf("I'm am Zombie\n"); break;}  // removes bad events where -nan is filled into the fMapped array
			fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1] = ii; 
		}

		for (Int_t ii=0; ii<nHits; ii++)
		{
        		if (ZombieAreAlive == kTRUE) {break;}

			counter_tot_tot++;
        		auto mapped = (R3BPspxMappedData const*)fMapped->At(ii);
        		//printf("Event Number: %d, Det: %d, Face: %d, Side: %d, Strip: %d, nHits: %d, Index: %d\n",counter_tot,mapped->GetDetector(),mapped->GetFace(),mapped->GetSide(),mapped->GetStrip(),nHits,fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1]);
    			if (MatchedEvents(ii) == kFALSE) {counter_matched_bad++; continue;}

			//Used to get the array index for the opposite side and opposite face of a detector
        		Int_t OtherSide = 0; Int_t OtherFace = 0;
			if (mapped->GetSide() == 1) {OtherSide = 2;}
			else if (mapped->GetSide() == 2) {OtherSide = 1;}
			if (mapped->GetFace() == 1) {OtherFace = 2;}
			else if (mapped->GetFace() == 2) {OtherFace = 1;}
			auto mapped2 = (R3BPspxMappedData const*)fMapped->At(fArray[mapped->GetDetector()-1][mapped->GetFace()-1][OtherSide-1][mapped->GetStrip()-1]);	

			Float_t time_corr1=0; Float_t time_corr2=0; Float_t time_corr_trig=0.;

			//need to add code so this is just for s473+s444 (Nik's FW)
        		if (fNumExpt == 473 || fNumExpt == 444)
			{
				time_corr1 = GetTimeCorr_s473(mapped->GetTime());
        			time_corr2 = GetTimeCorr_s473(mapped2->GetTime());
				time_corr_trig = GetTimeCorr_s473(mapped->GetTrigger());
			}

			//How to get the timing for s515 (CALIFA FW) 
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
					//printf("Det: %d, Face: %d, Side: %d, Strip:%d, Timing: %lf, Timing2: %lf\n",mapped->GetDetector()-1,mapped->GetFace()-1,mapped->GetSide()-1,mapped->GetStrip()-1,time_corr1,time_corr2);
					hTime[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1]->Fill(time_corr1);
					hEnergy[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-1]->Fill(mapped->GetEnergy());
					if (mapped->GetSide() == 1) hEnergyTot[(mapped->GetDetector()-1)][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mapped->GetEnergy() + mapped2->GetEnergy());
				}
				else {counter_mult_mismatch++;}

				hMult[mapped->GetDetector()-1]->Fill(MultFace(mapped->GetDetector()-1,mapped->GetFace()-1),MultFace(mapped->GetDetector()-1,OtherFace-1));

				//Interstrip Events - used for the position calibration for s473, s444 and s515
				if (mapped->GetSide()==1 && InterstripEvent(ii) != -1) 
				{
					Int_t interstrip_index;
					Bool_t WrongEnergy = kFALSE;
					//if (InterstripEvent(ii) >= 1000) {interstrip_index = InterstripEvent(ii)-1000;}
					//else {interstrip_index = InterstripEvent(ii);}
					interstrip_index = InterstripEvent(ii);
					
					//Get info info for each side of the neuighboring strip 
					auto mappedstrip1 = (R3BPspxMappedData const*)fMapped->At(interstrip_index);
					auto mappedstrip2 = (R3BPspxMappedData const*)fMapped->At(fArray[mappedstrip1->GetDetector()-1][mappedstrip1->GetFace()-1][OtherSide-1][mappedstrip1->GetStrip()-1]);
					//printf("InterstripIndex: %d, %d, partner: %d\n",InterstripEvent(ii),interstrip_index,fArray[mappedstrip1->GetDetector()-1][mappedstrip1->GetFace()-1][1][mappedstrip1->GetStrip()-1]);
					//printf("Mult: %d, Strip1 #%d:, Strip2 #:%d\n",MultFace(mapped->GetDetector()-1,mapped->GetFace()-1),mapped->GetStrip(),mappedstrip1->GetStrip());
				
					//if (mapped->GetDetector() == 5 || mapped->GetDetector() == 6) {printf("Interstrip. Detector: %d, index: %d\n",mapped->GetDetector(),InterstripEvent(ii));}
					
					//Total_Energy of 4 signals. Should equal approx. beam. If double beam, then it's 2 neighboring hits, not interstrip event
					Int_t totEnergy = mapped->GetEnergy() + mapped2->GetEnergy() + mappedstrip1->GetEnergy() + mappedstrip2->GetEnergy();
					Int_t Energy_strip1 = mapped->GetEnergy() + mapped2->GetEnergy();
					Int_t Energy_strip2 = mappedstrip1->GetEnergy() + mappedstrip2->GetEnergy();
						
					//if (InterstripEvent(ii) < 1000 && (totEnergy > maxenergy[mapped->GetDetector()-1] || totEnergy < maxenergy[mapped->GetDetector()-1]/100)) {WrongEnergy = kTRUE;}  
					if ((totEnergy > maxenergy[mapped->GetDetector()-1] || totEnergy < maxenergy[mapped->GetDetector()-1]/10)) {WrongEnergy = kTRUE; counter_bad_energy++;}  

					//Get info for corresponding back partner
					if (FacePartner(ii) != -1)
					{
						auto mapped3 = (R3BPspxMappedData const*)fMapped->At(FacePartner(ii));	

						//printf("FacePartner: %d, nHits: %d\n",FacePartner(ii),nHits);
						//Important to get calibration params for each strip on left and right side
						Double_t strip_id = -1, strip_back = -1;
						Bool_t IsRight = kFALSE; 
						if (mappedstrip1->GetStrip() > mapped->GetStrip()) {strip_id = (Double_t)mapped->GetStrip(); strip_back = (Double_t)mapped3->GetStrip(); IsRight = kTRUE;}
						else {strip_id = (Double_t)mapped->GetStrip() - 0.5; strip_back = (Double_t)mapped3->GetStrip() - 0.5; IsRight = kFALSE;}


						if (mapped->GetFace()==2) 
						{
							Double_t posx = ((Double_t)mapped->GetEnergy() - (Double_t)mapped2->GetEnergy()) / ((Double_t)mapped->GetEnergy() + (Double_t)mapped2->GetEnergy());
							Double_t posx1 = ((Double_t)mappedstrip1->GetEnergy() - (Double_t)mappedstrip2->GetEnergy()) / ((Double_t)mappedstrip1->GetEnergy() + (Double_t)mappedstrip2->GetEnergy());
							Double_t posy = (Double_t)mapped3->GetStrip();
							if ((posx-posx1) < -0.01 || (posx-posx1) > 0.01) {WrongEnergy = kTRUE;}

							if (InterstripEvent(ii) < 1000 && WrongEnergy == kFALSE)
							{
								counter_interstrip++; 
								hHitMap[(mapped->GetDetector()-1)]->Fill(mapped->GetStrip(),mapped3->GetStrip());
								//printf("Back Strip: %d, posx: %lf\n",mapped3->GetStrip(),posx);
								hPosX_det[(mapped->GetDetector()-1)]->Fill(posx,strip_back);
								hPosX[(mapped->GetDetector()-1)][(mapped->GetStrip()-1)]->Fill(strip_back,posx);
								hPosXvsStrip_det[(mapped->GetDetector()-1)]->Fill(strip_id,posx);
								hInterstrip[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mapped->GetEnergy());
								hInterstripTot[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(totEnergy);
								hStrip1vsStrip2[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(posx,posx1);
								if (IsRight == kTRUE) {hPosX_strip[mapped->GetDetector()-1][2*(mapped3->GetStrip()-1)]->Fill(posx);}
								else {hPosX_strip[mapped->GetDetector()-1][2*mapped3->GetStrip()-1]->Fill(posx);}
							}
						}
						else if (mapped->GetFace()==1) 
						{
							Double_t posx = (Double_t)mapped3->GetStrip();
							Double_t posy = ((Double_t)mapped2->GetEnergy() - (Double_t)mapped->GetEnergy()) / ((Double_t)mapped->GetEnergy() + (Double_t)mapped2->GetEnergy());
							Double_t posy1 = ((Double_t)mappedstrip2->GetEnergy() - (Double_t)mappedstrip1->GetEnergy()) / ((Double_t)mappedstrip1->GetEnergy() + (Double_t)mappedstrip2->GetEnergy());
							if ((posy-posy1) < -0.01 || (posy-posy1) > 0.01) {WrongEnergy = kTRUE;}
							
							if (InterstripEvent(ii) < 1000 && WrongEnergy == kFALSE)
							{
								counter_interstrip++; 
								hHitMap[(mapped->GetDetector()-1)]->Fill(mapped3->GetStrip(),mapped->GetStrip());
								//printf("Back Strip: %d, posy: %lf\n",mapped3->GetStrip(),posy);
								hPosY_det[(mapped->GetDetector()-1)]->Fill(strip_back,posy);
								hPosY[(mapped->GetDetector()-1)][(mapped->GetStrip()-1)]->Fill(strip_back,posy);
								hPosYvsStrip_det[(mapped->GetDetector()-1)]->Fill(strip_id,posy);
								hInterstrip[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(mapped->GetEnergy());
								hInterstripTot[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(totEnergy);
								hStrip1vsStrip2[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetStrip()-1]->Fill(posy,posy1);
								if (IsRight == kTRUE) {hPosY_strip[mapped->GetDetector()-1][2*(mapped3->GetStrip()-1)]->Fill(posy);}
								else {hPosY_strip[mapped->GetDetector()-1][2*mapped3->GetStrip()-1]->Fill(posy);}
							}
						}
					}
				}
       			}
       			else {counter_timing_bad++;}
		}
	}
	fMapped->Clear();
}


void R3BPspxMapped2CalPar::FinishEvent() 
{
	for (Int_t ii=0; ii<fNumDet; ii++)
	{
		for (Int_t jj=0; jj<2; jj++)
		{
			for (Int_t kk=0; kk<2; kk++)
			{
				for (Int_t ll=0; ll<fNumStrips; ll++)
				{
					fArray[ii][jj][kk][ll]=-1;
				}
			}
		}
	}
}

void R3BPspxMapped2CalPar::Reset() {}

void R3BPspxMapped2CalPar::GetPosParameters() 
{
	//do all the calibration parameter determination for the TGraphs here
	for (Int_t ii=0; ii<fNumDet; ii++)
	{
		for (Int_t jj=0; jj<fNumStrips; jj++)
		{
			if (hPosX[ii][jj]->GetEntries() > 0)
			{
				Int_t xcounter1=0, xcounter2=0;
				TProfile *profx1 = hPosX[ii][jj]->ProfileX();
				//TGraphErrors *graph1 = new TGraphErrors();
				//TGraphErrors *graph3 = new TGraphErrors();
				//TGraphErrors *graphres1 = new TGraphErrors();
				//TGraphErrors *graphres3 = new TGraphErrors();
				//printf("\n\nDetector: %d, Strip %d, Number of bins: %d\n",ii,jj,profx1->GetNbinsX());
			
				Double_t ratioX_left=0, ratioX_right=0, ratioX=0;
				
				for ( int ibin=1; ibin<profx1->GetNbinsX(); ibin++ )
				{
					double XX = profx1->GetBinCenter(ibin);
					if (profx1->GetBinEntries(ibin) > 20)
					{
						double YY = profx1->GetBinContent(ibin);
						double Error = profx1->GetBinError(ibin);

						double idealpos=0;
						double actualpos=0;
						double ratio=0;
						double residual=0;

						if (ibin % 2 == 0) 
						{
							idealpos = ((XX+0.25) - 16.)/16.;
							residual = YY - idealpos;
							actualpos = 16.*(YY+1.);
							ratio = actualpos/(XX+0.25);
							ratioX_right = ratio;

							//graph1->SetPoint(xcounter1,XX+0.25,YY); 
							//graph1->SetPointError(xcounter1,Error,Error); 
							//graphres1->SetPoint(xcounter1,XX+0.25,residual);

							xcounter1++;
						}
						else 
						{
							idealpos = ((XX-0.25) - 16.)/16.;
							residual = YY - idealpos;
							actualpos = 16.*(YY+1.);
							ratio = actualpos/(XX-0.25);
							ratioX_left = ratio;

							//graph3->SetPoint(xcounter2,XX-0.25,YY); 
							//graph3->SetPointError(xcounter2,Error,Error); 
							//graphres3->SetPoint(xcounter2,XX-0.25,residual); 

							xcounter2++;
						}
						if (ibin %2 == 0)
						{
							if (ratioX_left != 0 && ratioX_right != 0) {ratioX = (ratioX_left + ratioX_right)/2.;}
							else {ratioX = (ratioX_left + ratioX_right);}
							GainFactors[ii][0][jj][(Int_t)(XX+0.25)]=ratioX;
							outputfile << ii << "\t" << 0 << "\t" << jj << "\t" << (Int_t)(XX+0.25) << "\t" << ratioX << std::endl;
						}
						//printf("Inside Loop - XX: %lf, YY: %lf, idealpos: %lf, actualpos: %lf, residual: %lf, ratio: %lf, ii: %d\n",XX,YY,idealpos,actualpos,residual,XX/actualpos,ibin);
					}
					ratioX_left = 0; ratioX_right = 0; ratioX = 0;
				} // end loop over bins
			}

			
			if (hPosY[ii][jj]->GetEntries() > 0)
			{
				Int_t ycounter1=0, ycounter2=0;
				TProfile *profx2 = hPosY[ii][jj]->ProfileX();
				TGraphErrors *graph2 = new TGraphErrors();
				TGraphErrors *graph4 = new TGraphErrors();
				TGraphErrors *graphres2 = new TGraphErrors();
				TGraphErrors *graphres4 = new TGraphErrors();
				printf("\n\nDetector: %d, Strip %d, Number of bins: %d\n",ii,jj,profx2->GetNbinsX());
				
				Double_t ratioY_left=0, ratioY_right=0, ratioY=0;
				
				for ( int ibin=1; ibin<profx2->GetNbinsX(); ibin++ )
				{
					double XX = profx2->GetBinCenter(ibin);
					if (profx2->GetBinEntries(ibin) > 20)
					{
						double YY = profx2->GetBinContent(ibin);
						double Error = profx2->GetBinError(ibin);

						double idealpos=0;
						double actualpos=0;
						double ratio=0;
						double residual=0;

						if (ibin % 2 == 0) 
						{
							idealpos = ((XX+0.25) - 16.)/16.;
							residual = YY - idealpos;
							actualpos = 16.*(YY+1.);
							ratio = actualpos/(XX+0.25);
							ratioY_right = ratio;

							graph2->SetPoint(ycounter1,XX+0.25,YY); 
							graph2->SetPointError(ycounter1,Error,Error); 
							graphres2->SetPoint(ycounter1,XX+0.25,residual); 
						
							ycounter1++;
						}
						else 
						{
							idealpos = ((XX-0.25) - 16.)/16.;
							residual = YY - idealpos;
							actualpos = 16.*(YY+1.);
							ratio = actualpos/(XX-0.25);
							ratioY_left = ratio;

							graph4->SetPoint(ycounter2,XX-0.25,YY); 
							graph4->SetPointError(ycounter2,Error,Error); 
							graphres4->SetPoint(ycounter2,XX-0.25,residual); 
						
							ycounter2++;
						}
						if (ibin %2 == 0)
						{
							if (ratioY_right!=0 && ratioY_left!=0) {ratioY = (ratioY_right + ratioY_left)/2.;}
							else {ratioY = (ratioY_right + ratioY_left);}
							GainFactors[ii][1][jj][(Int_t)(XX+0.25)]=ratioY;
							outputfile << ii << "\t" << 1 << "\t" << jj << "\t" << (Int_t)(XX+0.25) << "\t" << ratioY << std::endl;
							//printf("Det: %d, Face: %d, Strip: %d, Backstrip: %lf, ratio: %lf\n",ii,1,jj,XX+0.25,ratioY);
						}
						//printf("Inside Loop - XX: %lf, YY: %lf, idealpos: %lf, actualpos: %lf, residual: %lf, ratio: %lf, ii: %d\n",XX,YY,idealpos,actualpos,residual,XX/actualpos,ibin);
					}
					ratioY_left = 0; ratioY_right = 0; ratioY = 0;
				} // end loop over bins
			}
		}
	}
}

void R3BPspxMapped2CalPar::FinishTask() 
{
	GetPosParameters();
	//for (Int_t ii=0; ii < fNumDet; ii++)
	for (Int_t ii=0; ii < 6; ii++)
	{
		Position->cd();
		hMult[ii]->Write();
		hPosX_det[ii]->Write();
		hPosY_det[ii]->Write();
		hHitMap[ii]->Write();
		hPosXvsStrip_det[ii]->Write();
		hPosYvsStrip_det[ii]->Write();

		Interstrip->cd();
		hTimeDiffPartner[ii]->Write();
		for (Int_t jj=0; jj < 2; jj++)
		{
			for (Int_t kk=0; kk < 2; kk++)
			{
				for (Int_t ll=0; ll < fNumStrips; ll++)
				{
					Timing->cd();
					hTime[ii][jj][kk][ll]->Write();
					Energy->cd();
					hEnergy[ii][jj][kk][ll]->Write();
					if (kk == 0)
					{
						hEnergyTot[ii][jj][ll]->Write();
						Interstrip->cd();
						hInterstrip[ii][jj][ll]->Write();
						hInterstripTot[ii][jj][ll]->Write();
						hStrip1vsStrip2[ii][jj][ll]->Write();
						
						PositionStrips->cd();
						hPosX_strip[ii][2*ll+jj]->Write();
						hPosY_strip[ii][2*ll+jj]->Write();

						if (jj == 0)
						{
							Position->cd();
							hPosX[ii][ll]->Write();
							hPosY[ii][ll]->Write();
						}
						//for (Int_t mm = 0; mm < fNumStrips; mm++)
					}
				}
			}
		}
	}
        printf("total events: %d, tot_tot events: %d, events of interest: %d, bad matched: %d, mulit mismatch: %d,bad timing: %d, interstrip: %d, bad_interstrip: %d, tot_interstrip: %d, bad_energy_partner: %d, bad_time_partner: %d, good_partner: %d, partner: %d, bad energy: %d\n",counter_tot, counter_tot_tot,counter_events, counter_matched_bad, counter_mult_mismatch, counter_timing_bad, counter_interstrip, counter_bad_interstrip,counter_tot_interstrip,counter_bad_time_partner,counter_bad_energy_partner,counter_good_partner,counter_partner,counter_bad_energy);
	outputfile.close();
}

Float_t R3BPspxMapped2CalPar::GetTimeCorr_s473(Int_t time)
{
	//corrects for the timing offset
	Int_t corr_time = time - 2048;
	if (corr_time < 0) {corr_time = 4096 + corr_time;}
	return corr_time;
}

Float_t R3BPspxMapped2CalPar::GetTimeCorr_s515(Int_t time, Int_t trigger)
{
	//corrects for the timing offset
	Int_t corr_time = time - trigger;
	return corr_time;
}

Bool_t R3BPspxMapped2CalPar::MatchedEvents(Int_t index) 
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
	else return kFALSE;
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

Int_t R3BPspxMapped2CalPar::FacePartner(Int_t index) 
{
	auto mapped = (R3BPspxMappedData const*)fMapped->At(index);

	Int_t Time=-1000; Int_t Time_partner=-9000;

	counter_partner++; 
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
			if (TotalEnergy > 0.2*maxenergy[mapped->GetDetector()-1]) {

				Time_partner = GetTimeCorr_s473(mappedpartner->GetTime()); 
				Time = GetTimeCorr_s473(mapped->GetTime());
				if ((Time - Time_partner) < 15 && (Time - Time_partner) > -15) 
				{
					faceindex = fArray[mapped->GetDetector()-1][OtherFace-1][mapped->GetSide()-1][jj]; 
					counter_good_partner++; 
					hTimeDiffPartner[mapped->GetDetector()-1]->Fill(GetTimeCorr_s473(mapped->GetTime())-GetTimeCorr_s473(mappedpartner->GetTime()));
				}
				else {counter_bad_time_partner++; continue;}
			}
			else {counter_bad_energy_partner++; continue;}
		}
	}
	return faceindex;
}

Int_t R3BPspxMapped2CalPar::MultFace(Int_t det, Int_t face)
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


//Find events with mult = 2 on 1 side and mult = 1 on other side. Make sure these events are neighbors. Make sure time difference is acceptable. Output is the index of the neighboring strip, same side 
Int_t R3BPspxMapped2CalPar::InterstripEvent(Int_t index) 
{
	auto mapped = (R3BPspxMappedData const*)fMapped->At(index);
	Int_t index_IS = -1;

	Int_t OtherFace=-1;
	if (mapped->GetFace()==1) OtherFace=2;
	else if (mapped->GetFace()==2) OtherFace=1;

	//if ((MultFace(mapped->GetDetector()-1, mapped->GetFace()-1)) != MultFace(mapped->GetDetector()-1,OtherFace-1)) {return -1;}
	if (MultFace(mapped->GetDetector()-1, mapped->GetFace()-1) != 2) {return -1;}
	else if (MultFace(mapped->GetDetector()-1, OtherFace-1) != 1) {return -1;}

	else 
	{
		//if (mapped->GetDetector()==5 || mapped->GetDetector()==6) {printf("Interstrip Event. Det: %d, Timing1: %d, TotEnergy: %d\n",mapped->GetDetector(),mapped->GetTime(),mapped->GetEnergy());}
		if (fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-2] == -1 && fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()] == -1) {return -1;}
		//else if (fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-2] != -1 && fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()] != -1) {return -1;}
		else if (fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-2] != -1) {index_IS = fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()-2];}
		else if (fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()] != -1)   {index_IS = fArray[mapped->GetDetector()-1][mapped->GetFace()-1][mapped->GetSide()-1][mapped->GetStrip()];}
		else return -1;

		if (index_IS < 0 || index_IS > fMapped->GetEntries()) {return -1;}

		if (index_IS != -1 && MatchedEvents(index_IS) == kTRUE) 
		{
			counter_tot_interstrip++;
			auto mapped2 = (R3BPspxMappedData const*)fMapped->At(index_IS);
			//if (mapped->GetDetector()==5 || mapped->GetDetector()==6) {printf("Interstrip Event. Det: %d, Timing1: %d, Timing2: %d, TotEnergy: %d\n",mapped->GetDetector(),mapped->GetTime(),mapped2->GetTime(),(mapped->GetEnergy()+mapped2->GetEnergy()));}
			//if ((mapped->GetTime()-mapped2->GetTime()) < 15 && (mapped->GetTime()-mapped2->GetTime()) > -15) 
			if ((mapped->GetTime()-mapped2->GetTime()) < 5 && (mapped->GetTime()-mapped2->GetTime()) > -5) 
			{
				return index_IS; 
			}
			else 
			{
				counter_bad_interstrip++;
				return -1;
				//return index_IS+1000;
			}
		}
	}
	return -1;
}

void R3BPspxMapped2CalPar::SetNumDet(Int_t det) {fNumDet = det;}
void R3BPspxMapped2CalPar::SetNumStrips(Int_t strip) {fNumStrips = strip;}
void R3BPspxMapped2CalPar::SetNumExpt(Int_t expt) {fNumExpt = expt;}
void R3BPspxMapped2CalPar::SetParOutName(TString paroutname) {fParOutName = paroutname;}

ClassImp(R3BPspxMapped2CalPar)
