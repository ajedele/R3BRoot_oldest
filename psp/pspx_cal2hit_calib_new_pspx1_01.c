/*
   An attempt at a new script for the last calibration step, cal2hit, for the PSP detectors.
   written by Ashton Falduto 27.07.2021
 */
#include <TROOT.h>
#include "TStyle.h"
#include "TMath.h"
#include <TChain.h>
#include <TFile.h>
#include "TObject.h"
#include "TClonesArray.h"

#include "TCanvas.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TInterpreter.h"

#include "TF1.h"
#include "TText.h"
#include "FairEventHeader.h"
#include "/u/ekudaib/R3BRootnew/R3BRoot/r3bbase/R3BEventHeader.h"
#include <iostream>

using namespace std;
using namespace TMath;
void pspx_cal2hit_test(UInt_t run_min = 508, UInt_t run_max = 508, Double_t threshold=0, bool save=false,bool check = false)
{
	UInt_t verbosity_step_size = 100000; // number, how often progress of skript is reported

	// Detector specific variables
	UInt_t nStrips = 32;
	UInt_t npeaks = nStrips - 1;
	Double_t strip_thickness = 9.57 / nStrips; //strip thickness was adapted to R3BWiki value
	Double_t rangeE = 2000; // cut range for same energy on front and back, was changed from 18500
	Double_t neve=1;

	TChain *fChain = new TChain("evt");
	TString filename;
	for(UInt_t i=run_min; i<=run_max; i++)
	{
		filename = Form("/lustre/land/ekudaib/rootfiles/s515/main_run%03d_r3broot_precal2cal.root", i);
		cout << "Open " << filename << endl;
		fChain->Add(filename);
	}
	#include "new_class.h"
	fChain->SetMakeClass(1);

	//Calibration parameters
	Double_t mean_strip_X[nStrips];
	Double_t mean_strip_Y[nStrips];
	Double_t pos_strip[nStrips];
	Double_t slope_strip_X[nStrips];
	Double_t offset_strip_X[nStrips];
	Double_t slope_strip_Y[nStrips];
	Double_t offset_strip_Y[nStrips];

	for(UInt_t i=0;i<=nStrips;i++)
	{
		mean_strip_X[i]=0;
		mean_strip_Y[i]=0;
		pos_strip[i]=(i*strip_thickness-nStrips/2.*strip_thickness);
		slope_strip_X[i]=0;
		offset_strip_X[i]=strip_thickness*nStrips/2;
		slope_strip_Y[i]=0;
		offset_strip_Y[i]=strip_thickness*nStrips/2;
		cout<<i+1<<" pos: "<<pos_strip[i]<<endl;
	}

	// --- HISTOGRAMS ---
	Double_t binPos = 500;

	Double_t minPos = -0.5;
	Double_t maxPos = 0.5;

	TH1F* hPosX;
	TH1F* hPosY;
	TGraph* hCalibX;
	TGraph* hCalibY;
			
	hPosX = new TH1F("H_pos_strip_X","Position X, X",binPos,minPos,maxPos);
	hPosX->GetXaxis()->SetTitle("Position/a.u.");
	hPosX->GetYaxis()->SetTitle("Counts");

	hPosY= new TH1F(Form("pos_Y"),Form("Position Y, Y"),binPos,minPos,maxPos);
	hPosY->GetXaxis()->SetTitle("Position/a.u.");
	hPosY->GetYaxis()->SetTitle("Counts");

	// --- EVENTLOOP ---
	if (fChain == 0)
	return;
	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;

	Int_t multx, multy;
	Double_t E_X, E_Y;
	Int_t strip_X[nStrips];
	Int_t strip_Y[nStrips];

	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Pspx1_xCal*",1);
	fChain->SetBranchStatus("Pspx1_yCal*",1);

	for (Long64_t jentry=0; jentry<neve*nentries;jentry++)
	{
		Long64_t ientry = fChain->LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		if(nb==4) continue;
		if(nb==0) {cout<< "Something went terribly wrong!" << endl; break;};
		if (jentry%verbosity_step_size==0)
		{
			cout << jentry << " of " << nentries << " == "<< 100. * jentry / nentries <<"% of the events analyzed\r"<< flush;
		}

		//Reseting all values
		multx = 0;
		multy = 0;
		E_X = 0;
		E_Y = 0;

		for(UInt_t i=0;i<nStrips;i++)
		{
			strip_X[i]=0;
			strip_Y[i]=0;
		}

		//Get the multiplicty for every side of the detector
		for(Int_t i=0; i<Pspx1_xCal_;i++)
		{
			if(abs(Pspx1_xCal_fEnergy[i])>threshold)
			{
				strip_X[multx]=Pspx1_xCal_fStrip[i];
				multx++;
				E_X += Pspx1_xCal_fEnergy[i];
			}
		}
		for(Int_t i=0; i<Pspx1_yCal_;i++)
		{
			if(abs(Pspx1_yCal_fEnergy[i])>threshold)
			{
				strip_Y[multy]=Pspx1_yCal_fStrip[i];
				E_Y += Pspx1_yCal_fEnergy[i];
				multy++;
			}
		}

		// middle interstrip event
		if (multx == 1 && multy == 2 && ((abs(E_X) < abs(E_Y - rangeE) && abs(E_X) > abs(E_Y + rangeE)) ||(abs(E_Y) < abs(E_X + rangeE) && abs(E_Y) > abs(E_X - rangeE))))
		{
			for(Int_t i=0; i<Pspx1_xCal_;i++)//pod voprosom
			{
				hPosX->Fill(Pspx1_xCal_fPos[i]);
			}
								
		}
		if (multx == 2 && multy == 1 && ((abs(E_X) < abs(E_Y - rangeE) && abs(E_X) > abs(E_Y + rangeE)) || (abs(E_Y) < abs(E_X + rangeE) && abs(E_Y) > abs(E_X - rangeE))))
		{
			for(Int_t i=0; i<Pspx1_yCal_;i++)
			{
				hPosY->Fill(Pspx1_yCal_fPos[i]);
			}
		}
	}
												
		//Calibration/ Fit

		TF1 * gaus_pol1 = new TF1("gaus_pol1","gaus(0)",-1,1);
		TF1 * gaus_pol2 = new TF1("gaus_pol2","gaus(0)",-1,1);
		TF1 * gaus_pol3 = new TF1("gaus_pol3","gaus(0)",-1,1);
		TF1 * gaus_pol4 = new TF1("gaus_pol4","gaus(0)",-1,1);
		TF1 * gaus_pol5 = new TF1("gaus_pol5","gaus(0)",-1,1);
		TF1 * gaus_pol6 = new TF1("gaus_pol6","gaus(0)",-1,1);
		TF1 * gaus_pol7 = new TF1("gaus_pol7","gaus(0)",-1,1);
		TF1 * gaus_pol8 = new TF1("gaus_pol8","gaus(0)",-1,1);
		/*TF1 * gaus_pol9 = new TF1("gaus_pol9","gaus(0)",-1,1);
		TF1 * gaus_pol10 = new TF1("gaus_pol10","gaus(0)",-1,1);
		TF1 * gaus_pol11 = new TF1("gaus_pol11","gaus(0)",-1,1);
		TF1 * gaus_pol12 = new TF1("gaus_pol12","gaus(0)",-1,1);
		TF1 * gaus_pol13 = new TF1("gaus_pol13","gaus(0)",-1,1);
		TF1 * gaus_pol14 = new TF1("gaus_pol14","gaus(0)",-1,1);
		TF1 * gaus_pol15 = new TF1("gaus_pol15","gaus(0)",-1,1);*/

		//TF1 * gaus_pol16 = new TF1("gaus_pol16","gaus(0)",-1,1);
		TF1 * gaus_pol17 = new TF1("gaus_pol17","gaus(0)",-1,1);
		TF1 * gaus_pol18 = new TF1("gaus_pol18","gaus(0)",-1,1);
		TF1 * gaus_pol19 = new TF1("gaus_pol19","gaus(0)",-1,1);
		TF1 * gaus_pol20 = new TF1("gaus_pol20","gaus(0)",-1,1);
		TF1 * gaus_pol21 = new TF1("gaus_pol21","gaus(0)",-1,1);
		TF1 * gaus_pol22 = new TF1("gaus_pol22","gaus(0)",-1,1);
		TF1 * gaus_pol23 = new TF1("gaus_pol23","gaus(0)",-1,1);
		TF1 * gaus_pol24 = new TF1("gaus_pol24","gaus(0)",-1,1);
		TF1 * gaus_pol25 = new TF1("gaus_pol25","gaus(0)",-1,1);
		TF1 * gaus_pol26 = new TF1("gaus_pol26","gaus(0)",-1,1);
		TF1 * gaus_pol27 = new TF1("gaus_pol27","gaus(0)",-1,1);
		TF1 * gaus_pol28 = new TF1("gaus_pol28","gaus(0)",-1,1);
		TF1 * gaus_pol29 = new TF1("gaus_pol29","gaus(0)",-1,1);
		//TF1 * gaus_pol30 = new TF1("gaus_pol30","gaus(0)",-1,1);	
	
		//X

		gaus_pol1->SetParLimits(0,1,1600);
		gaus_pol1->SetParLimits(2,0.001,0.25);
		gaus_pol1->SetParLimits(3,0,1000);

		gaus_pol2->SetParLimits(0,1,4500);
		gaus_pol2->SetParLimits(2,0.001,0.3);
		gaus_pol2->SetParLimits(3,0,1000);

		gaus_pol3->SetParLimits(0,1,4800);
		gaus_pol3->SetParLimits(2,0.0001,0.25);
		gaus_pol3->SetParLimits(3,0,1000);

		gaus_pol4->SetParLimits(0,1,5400);
		gaus_pol4->SetParLimits(2,0.00001,0.05);
		gaus_pol4->SetParLimits(3,0,1000);

		gaus_pol5->SetParLimits(0,1,5400);
		gaus_pol5->SetParLimits(2,0.001,0.2);
		gaus_pol5->SetParLimits(3,0,1000);

		gaus_pol6->SetParLimits(0,1,4600);
		gaus_pol6->SetParLimits(2,0.0001,0.1);
		gaus_pol6->SetParLimits(3,0,1000);

		gaus_pol7->SetParLimits(0,1,3000);
		gaus_pol7->SetParLimits(2,0.001,0.25);
		gaus_pol7->SetParLimits(3,0,1000);

		gaus_pol8->SetParLimits(0,1,500);
		gaus_pol8->SetParLimits(2,0.001,0.3);
		gaus_pol8->SetParLimits(3,0,1000);

		/*gaus_pol9->SetParLimits(0,1,5000);
		gaus_pol9->SetParLimits(2,0.0001,0.25);
		gaus_pol9->SetParLimits(3,0,1000);

		gaus_pol10->SetParLimits(0,1,5000);
		gaus_pol10->SetParLimits(2,0.00001,0.05);
		gaus_pol10->SetParLimits(3,0,1000);

		gaus_pol11->SetParLimits(0,1,5000);
		gaus_pol11->SetParLimits(2,0.00001,0.05);
		gaus_pol11->SetParLimits(3,0,1000);

		gaus_pol12->SetParLimits(0,1,5000);
		gaus_pol12->SetParLimits(2,0.00001,0.05);
		gaus_pol12->SetParLimits(3,0,1000);

		gaus_pol13->SetParLimits(0,1,5000);
		gaus_pol13->SetParLimits(2,0.001,0.2);
		gaus_pol13->SetParLimits(3,0,1000);

		gaus_pol14->SetParLimits(0,1,5000);
		gaus_pol14->SetParLimits(2,0.0001,0.1);
		gaus_pol14->SetParLimits(3,0,1000);


		gaus_pol15->SetParLimits(0,1,5000);
		gaus_pol15->SetParLimits(2,0.0001,0.1);
		gaus_pol15->SetParLimits(3,0,1000);*/

		//Y

		//gaus_pol16->SetParLimits(0,1,20);
		//gaus_pol16->SetParLimits(2,0.0001,0.1);
		//gaus_pol16->SetParLimits(3,0,1000);

		gaus_pol17->SetParLimits(0,1,180);
		gaus_pol17->SetParLimits(2,0.001,0.1);
		gaus_pol17->SetParLimits(3,0,1000);

		gaus_pol18->SetParLimits(0,1,450);
		gaus_pol18->SetParLimits(2,0.0001,0.35);
		gaus_pol18->SetParLimits(3,0,1000);

		gaus_pol19->SetParLimits(0,1,1000);
		gaus_pol19->SetParLimits(2,0.0001,0.25);
		gaus_pol19->SetParLimits(3,0,1000);

		gaus_pol20->SetParLimits(0,1,1420);
		gaus_pol20->SetParLimits(2,0.0001,0.25);
		gaus_pol20->SetParLimits(3,0,1000);

		gaus_pol21->SetParLimits(0,1,1220);
		gaus_pol21->SetParLimits(2,0.0001,0.25);
		gaus_pol21->SetParLimits(3,0,1000);

		gaus_pol22->SetParLimits(0,1,1610);
		gaus_pol22->SetParLimits(2,0.0001,0.1);
		gaus_pol22->SetParLimits(3,0,1000);

		gaus_pol23->SetParLimits(0,1,1750);
		gaus_pol23->SetParLimits(2,0.001,0.1);
		gaus_pol23->SetParLimits(3,0,1000);

		gaus_pol24->SetParLimits(0,1,1405);
		gaus_pol24->SetParLimits(2,0.0001,0.35);
		gaus_pol24->SetParLimits(3,0,1000);

		gaus_pol25->SetParLimits(0,1,1400);
		gaus_pol25->SetParLimits(2,0.0001,0.25);
		gaus_pol25->SetParLimits(3,0,1000);

		gaus_pol26->SetParLimits(0,1,1210);
		gaus_pol26->SetParLimits(2,0.0001,0.25);
		gaus_pol26->SetParLimits(3,0,1000);

		gaus_pol27->SetParLimits(0,1,960);
		gaus_pol27->SetParLimits(2,0.0001,0.25);
		gaus_pol27->SetParLimits(3,0,1000);
		
		gaus_pol28->SetParLimits(0,1,600);
		gaus_pol28->SetParLimits(2,0.0001,0.25);
		gaus_pol28->SetParLimits(3,0,1000);

		gaus_pol29->SetParLimits(0,1,350);
		gaus_pol29->SetParLimits(2,0.0001,0.25);
		gaus_pol29->SetParLimits(3,0,1000);

		//gaus_pol30->SetParLimits(0,1,200);
		//gaus_pol30->SetParLimits(2,0.0001,0.25);
		//gaus_pol30->SetParLimits(3,0,1000);

		TSpectrum *s = new TSpectrum(15);
		Int_t nfound = 0;

		Double_t *xpeaks;
		Double_t fit_range = 0.1;
		Double_t fit_min = -fit_range;
		Double_t fit_max = fit_range;

		nfound = s->Search(hPosX,2,"",0.05);
		xpeaks = s->GetPositionX();

		gaus_pol1->SetParameter(0,5000);
		gaus_pol1->SetParameter(1,xpeaks[0]);
		gaus_pol1->SetParameter(2,0.05);

		gaus_pol2->SetParameter(0,5000);
		gaus_pol2->SetParameter(1,xpeaks[0]);
		gaus_pol2->SetParameter(2,0.05);

		gaus_pol3->SetParameter(0,5000);
		gaus_pol3->SetParameter(1,xpeaks[0]);
		gaus_pol3->SetParameter(2,0.05);

		gaus_pol4->SetParameter(0,5000);
		gaus_pol4->SetParameter(1,xpeaks[0]);
		gaus_pol4->SetParameter(2,0.05);

		gaus_pol5->SetParameter(0,5000);
		gaus_pol5->SetParameter(1,xpeaks[0]);
		gaus_pol5->SetParameter(2,0.05);

		gaus_pol6->SetParameter(0,5000);
		gaus_pol6->SetParameter(1,xpeaks[0]);
		gaus_pol6->SetParameter(2,0.05);

		gaus_pol7->SetParameter(0,5000);
		gaus_pol7->SetParameter(1,xpeaks[0]);
		gaus_pol7->SetParameter(2,0.05);

		gaus_pol8->SetParameter(0,5000);
		gaus_pol8->SetParameter(1,xpeaks[0]);
		gaus_pol8->SetParameter(2,0.05);

		/*gaus_pol9->SetParameter(0,5000);
		gaus_pol9->SetParameter(1,xpeaks[0]);
		gaus_pol9->SetParameter(2,0.05);

		gaus_pol10->SetParameter(0,5000);
		gaus_pol10->SetParameter(1,xpeaks[0]);
		gaus_pol10->SetParameter(2,0.05);

		gaus_pol11->SetParameter(0,5000);
		gaus_pol11->SetParameter(1,xpeaks[0]);
		gaus_pol11->SetParameter(2,0.05);

		gaus_pol12->SetParameter(0,5000);
		gaus_pol12->SetParameter(1,xpeaks[0]);
		gaus_pol12->SetParameter(2,0.05);

		gaus_pol13->SetParameter(0,5000);
		gaus_pol13->SetParameter(1,xpeaks[0]);
		gaus_pol13->SetParameter(2,0.05);

		gaus_pol14->SetParameter(0,5000);
		gaus_pol14->SetParameter(1,xpeaks[0]);
		gaus_pol14->SetParameter(2,0.05);

		gaus_pol15->SetParameter(0,5000);
		gaus_pol15->SetParameter(1,xpeaks[0]);
		gaus_pol15->SetParameter(2,0.05);*/

		hPosX->Fit(gaus_pol1,"BIQ+","",-0.13,-0.11); //6set
		hPosX->Fit(gaus_pol2,"BIQ+","",-0.075,-0.045); //6set     
		hPosX->Fit(gaus_pol3,"BIQ+","",-0.02,0.02);  //6set   
		hPosX->Fit(gaus_pol4,"BIQ+","", 0.045,0.075); //6set   
		hPosX->Fit(gaus_pol5,"BIQ+","", 0.1,0.14);//6set  
		hPosX->Fit(gaus_pol6,"BIQ+","", 0.17,0.19);    //6set
		hPosX->Fit(gaus_pol7,"BIQ+","", 0.23,0.255); //6set
		hPosX->Fit(gaus_pol8,"BIQ+","", 0.29,0.31); //6set     
		/*hPosX->Fit(gaus_pol9,"BIQ+","",-0.075,-0.045);  //6set   
		hPosX->Fit(gaus_pol10,"BIQ+","",-0.015,0.015); //6set   
		hPosX->Fit(gaus_pol11,"BIQ+","",0.045,0.075);//6set  
		hPosX->Fit(gaus_pol12,"BIQ+","",0.105,0.1375);    //6set
		hPosX->Fit(gaus_pol13,"BIQ+","",-0.015,0.015); //6set   
		hPosX->Fit(gaus_pol14,"BIQ+","",0.045,0.075);//6set  
		hPosX->Fit(gaus_pol15,"BIQ+","",0.105,0.1375);    //6set*/

		mean_strip_X[0]=gaus_pol1->GetParameter(1);
		mean_strip_X[1]=gaus_pol2->GetParameter(1);
		mean_strip_X[2]=gaus_pol3->GetParameter(1);
		mean_strip_X[3]=gaus_pol4->GetParameter(1);
		mean_strip_X[4]=gaus_pol5->GetParameter(1);
		mean_strip_X[5]=gaus_pol6->GetParameter(1);
		mean_strip_X[6]=gaus_pol7->GetParameter(1);
		mean_strip_X[7]=gaus_pol8->GetParameter(1);
		/*mean_strip_X[8]=gaus_pol9->GetParameter(1);
		mean_strip_X[9]=gaus_pol10->GetParameter(1);
		mean_strip_X[10]=gaus_pol11->GetParameter(1);
		mean_strip_X[11]=gaus_pol12->GetParameter(1);
		mean_strip_X[12]=gaus_pol13->GetParameter(1);
		mean_strip_X[13]=gaus_pol14->GetParameter(1);
		mean_strip_X[14]=gaus_pol15->GetParameter(1);*/
		
	
		nfound = s->Search(hPosY,2,"",0.05);
		xpeaks = s->GetPositionX();


		//gaus_pol16->SetParameter(0,5000);
		//gaus_pol16->SetParameter(1,xpeaks[0]);
		//gaus_pol16->SetParameter(2,0.05);

		gaus_pol17->SetParameter(0,5000);
		gaus_pol17->SetParameter(1,xpeaks[0]);
		gaus_pol17->SetParameter(2,0.05);

		gaus_pol18->SetParameter(0,5000);
		gaus_pol18->SetParameter(1,xpeaks[0]);
		gaus_pol18->SetParameter(2,0.05);

		gaus_pol19->SetParameter(0,5000);
		gaus_pol19->SetParameter(1,xpeaks[0]);
		gaus_pol19->SetParameter(2,0.05);

		gaus_pol20->SetParameter(0,5000);
		gaus_pol20->SetParameter(1,xpeaks[0]);
		gaus_pol20->SetParameter(2,0.05);

		gaus_pol21->SetParameter(0,5000);
		gaus_pol21->SetParameter(1,xpeaks[0]);
		gaus_pol21->SetParameter(2,0.05);

		gaus_pol22->SetParameter(0,5000);
		gaus_pol22->SetParameter(1,xpeaks[0]);
		gaus_pol22->SetParameter(2,0.05);

		gaus_pol23->SetParameter(0,5000);
		gaus_pol23->SetParameter(1,xpeaks[0]);
		gaus_pol23->SetParameter(2,0.05);

		gaus_pol24->SetParameter(0,5000);
		gaus_pol24->SetParameter(1,xpeaks[0]);
		gaus_pol24->SetParameter(2,0.05);

		gaus_pol25->SetParameter(0,5000);
		gaus_pol25->SetParameter(1,xpeaks[0]);
		gaus_pol25->SetParameter(2,0.05);

		gaus_pol26->SetParameter(0,5000);
		gaus_pol26->SetParameter(1,xpeaks[0]);
		gaus_pol26->SetParameter(2,0.05);

		gaus_pol27->SetParameter(0,5000);
		gaus_pol27->SetParameter(1,xpeaks[0]);
		gaus_pol27->SetParameter(2,0.05);

		gaus_pol28->SetParameter(0,5000);
		gaus_pol28->SetParameter(1,xpeaks[0]);
		gaus_pol28->SetParameter(2,0.05);

		gaus_pol29->SetParameter(0,5000);
		gaus_pol29->SetParameter(1,xpeaks[0]);
		gaus_pol29->SetParameter(2,0.05);

		//gaus_pol30->SetParameter(0,5000);
		//gaus_pol30->SetParameter(1,xpeaks[0]);
		//gaus_pol30->SetParameter(2,0.05);

		//hPosY->Fit(gaus_pol16,"BIQ+","",-0.44,-0.42); //6set
		hPosY->Fit(gaus_pol17,"BIQ+","",-0.386,-0.36); //6set     
		hPosY->Fit(gaus_pol18,"BIQ+","",-0.328,-0.3);  //6set   
		hPosY->Fit(gaus_pol19,"BIQ+","",-0.252,-0.238); //6set   
		hPosY->Fit(gaus_pol20,"BIQ+","",-0.19,-0.178);//6set  
		hPosY->Fit(gaus_pol21,"BIQ+","",-0.13,-0.114);    //6set
		hPosY->Fit(gaus_pol22,"BIQ+","",-0.07,-0.05); //6set
		hPosY->Fit(gaus_pol23,"BIQ+","",-0.016,0.008); //6set   

		hPosY->Fit(gaus_pol24,"BIQ+","",0.055,0.075);  //6set  
		hPosY->Fit(gaus_pol25,"BIQ+","",0.11,0.1325);  //6set 
		hPosY->Fit(gaus_pol26,"BIQ+","",0.18,0.195); //6set   
		hPosY->Fit(gaus_pol27,"BIQ+","",0.24,0.26);//6set  
		hPosY->Fit(gaus_pol28,"BIQ+","",0.3,0.325);    //6set
		hPosY->Fit(gaus_pol29,"BIQ+","",0.362,0.375); //6set   
		//hPosY->Fit(gaus_pol30,"BIQ+","",0.42,0.44);//6set  

		//mean_strip_Y[0]=gaus_pol16->GetParameter(1);
		mean_strip_Y[0]=gaus_pol17->GetParameter(1);
		mean_strip_Y[1]=gaus_pol18->GetParameter(1);
		mean_strip_Y[2]=gaus_pol19->GetParameter(1);
		mean_strip_Y[3]=gaus_pol20->GetParameter(1);
		mean_strip_Y[4]=gaus_pol21->GetParameter(1);
		mean_strip_Y[5]=gaus_pol22->GetParameter(1);
		mean_strip_Y[6]=gaus_pol23->GetParameter(1);
		mean_strip_Y[7]=gaus_pol24->GetParameter(1);
		mean_strip_Y[8]=gaus_pol25->GetParameter(1);
		mean_strip_Y[9]=gaus_pol26->GetParameter(1);
		mean_strip_Y[10]=gaus_pol27->GetParameter(1);
		mean_strip_Y[11]=gaus_pol28->GetParameter(1);
		mean_strip_Y[12]=gaus_pol29->GetParameter(1);
		//mean_strip_Y[14]=gaus_pol30->GetParameter(1);
		

		TF1 * lin_pol1 = new TF1("lin_pol1","pol1(0)",-nStrips/2.*strip_thickness, nStrips/2.*strip_thickness);
		Double_t pos_X[8];
		Double_t pos_calib_X[8];
		Double_t pos_Y[13];
		Double_t pos_calib_Y[13];
							
		for(Int_t i=0;i<8;i++)
		{
			pos_calib_X[i] = pos_strip[i+14];
			pos_X[i] = mean_strip_X[i];	
			cout<<i+1<<", "<<pos_X[i]<<", "<<pos_X[i]*4.785<<", "<<pos_calib_X[i]<<endl;						
        }
		
		for(Int_t i=0;i<13;i++)
		{
			pos_Y[i] = mean_strip_Y[i];	
			pos_calib_Y[i] = pos_strip[i+10];	
			cout<<i+1<<", "<<pos_Y[i]<<", "<<pos_Y[i]*4.785<<", "<<pos_calib_Y[i]<<endl;						
        }

							
        hCalibX=new TGraph(8,pos_X, pos_calib_X);
        hCalibX->GetXaxis()->SetTitle("Position/a.u.");
        hCalibX->GetYaxis()->SetTitle("Position/cm");
		hCalibX->Fit("lin_pol1","IQ+","",-1,1);

		offset_strip_X[0] = lin_pol1->GetParameter(0);
		slope_strip_X[0] = lin_pol1->GetParameter(1);

		hCalibY=new TGraph(13,pos_Y,pos_calib_Y);
        hCalibY->GetXaxis()->SetTitle("Position/a.u.");
        hCalibY->GetYaxis()->SetTitle("Position/cm");
		hCalibY->Fit("lin_pol1","IQ+","",-1,1);

		offset_strip_Y[0] = lin_pol1->GetParameter(0);
		slope_strip_Y[0] = lin_pol1->GetParameter(1);

		cout<<"offset front:"<<offset_strip_X[0]<<", "<<"slope front:"<<slope_strip_X[0]<<endl;
		cout<<"offset back:"<<offset_strip_Y[0]<<", "<<"slope back:"<<slope_strip_Y[0]<<endl;
							
		//Plots
							
		gStyle->SetPalette(1);
		gStyle->SetOptStat(110011);
		gStyle->SetOptFit(0011);
		gStyle->SetPalette(1);
				
		hPosX->SetLineColor(4);
		hPosY->SetLineColor(4);
							
		TCanvas *c1 = new TCanvas("PSPX1_X_center","PSPX1 X Center ",1400,700);
		c1->Divide(1,1);
		c1->cd(1);
		hPosX->Draw();

		TCanvas * c2 = new TCanvas("PSPX1_Y_center","PSPX1 Y Center ",1400,700);
		c2->Divide(1,1);
		c2->cd(1);
		hPosY->Draw();
							
		TCanvas * c3 = new TCanvas("TGraph","TGraph",1400,700);
		c3->Divide(1,2);
		c3->cd(1);
		hCalibX->Draw("AP*");
		c3->cd(2);
		hCalibY->Draw("AP*");
}						
		