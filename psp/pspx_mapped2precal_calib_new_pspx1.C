/*
 * Automated calibration for PSPX1 detectors.
 *
 * HOWTO use: (not yet)
 * run script to get Mapped2Precal calibration parameters ("root -l pspx_mapped2precal_calib_frontback.C+(true)")
 * and copy them into the corresponding parameter file
 *
 * WARNING: The arrays in this script start counting with 0 (i.e. 0-4 for detectors, 0-15 for strips, 0-64 for channels), but
 * the printed output is for an array gain starting with 1
 */

#include <TROOT.h>
#include "TApplication.h"
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

#include "TF1.h"

#include "/u/ajedele/R3BRoot_oldest/r3bbase/R3BEventHeader.h"
#include "FairEventHeader.h"

#include <iostream>

using namespace std;
using namespace TMath;

void pspx_mapped2precal_calib_new_pspx1(UInt_t run_min = 230, UInt_t run_max=230, bool save = false, bool check = false) {

	UInt_t verbosity_step_size = 100000; // number, how often progress of skript is reported

	// Detector specific variables
	UInt_t nStrips = 32;
	UInt_t npeaks = nStrips - 1;
	UInt_t nChannels = 4 * nStrips;

	Double_t strip_thickness = 9.57 / nStrips; //strip thickness was adapted to R3BWiki value

	Int_t strip_f0=17, strip_f1=18, strip_b0=16, strip_b1=17;
	Double_t rangeE = 10000000000; // cut range for same energy on front and back, was changed from 18500
	Double_t neve=1;
	TChain *fChain = new TChain("evt");
	TString filename;
	for(UInt_t i=run_min; i<=run_max; i++){
		if (check) {
			filename = Form("/home/matthias/Work/r3b/pspx/pspx_run%i_hit_hitpar.root", i);
		} else {
			// filename = Form("/home/matthias/Work/r3b/pspx/pspx_run%04d_hit_hitpar.root", i);
			filename = Form("/lustre/land/ajedele/s473/rootfiles/main%04d_pspx.root", i);
		}
		cout << "Open " << filename << endl;
		fChain->Add(filename);
	}

// #include "pspx_new.h"
#include "new_class.h"

	fChain->SetMakeClass(1);
	// Calibration parameters
	Double_t mean_strip[nStrips];
	Double_t gain_strip[nStrips];

	for (UInt_t i = 0; i < nStrips; i++) {
		mean_strip[i] = 0;
		gain_strip[i] = 0;
	}
	// --- HISTOGRAMS ---
	Double_t binPos = 500;
	
	Double_t minPos = -0.5;
	Double_t maxPos = 0.5;
	
	Double_t minPosFit = -0.1;
	Double_t maxPosFit = 0.1;

	Double_t titlesize = 0.045;

	TH1F* hPosX[nStrips];
	TH1F* hPosY[nStrips];

	TH1F* hTimeX[nStrips];
	TH1F* hTimeY[nStrips];
	TH1F* hTimeDiffX[nStrips];
	TH1F* hTimeDiffY[nStrips];

	for (UInt_t k = 0; k < nStrips; k++) {
		hPosX[k] =
			new TH1F(Form("pos_strip_x_middle%d", k), Form("Position, Strip %d, X; Position/arb.u.; Counts", k+1), binPos, minPosFit, maxPosFit);
		hPosY[k] =
			new TH1F(Form("pos_strip_y_middle%d", k), Form("Position, Strip %d, Y; Position/arb.u.; Counts", k+1) , binPos, minPosFit, maxPosFit);

		hPosX[k]->GetXaxis()->SetTitleSize(titlesize);
		hPosX[k]->GetYaxis()->SetTitleSize(titlesize);
		hPosX[k]->GetYaxis()->SetTitleOffset(1);

		hPosY[k]->GetXaxis()->SetTitleSize(titlesize);
		hPosY[k]->GetYaxis()->SetTitleSize(titlesize);
		hPosY[k]->GetYaxis()->SetTitleOffset(1);
	}
	TH2F *hinterstrip = new TH2F("hinterstrip","hinterstrip", binPos, minPos, maxPos, binPos, minPos, maxPos);
	TH2F *hstrip = new TH2F("hstrip","hstrip", 32,1,33,32,1,33);

	// --- EVENTLOOP ---
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;

	Int_t mult_front, mult_back;
	Double_t energy_front, energy_back;

	Int_t counter_front = 0; Int_t counter_back = 0;

	Int_t strip_front[nStrips]; Int_t strip_back[nStrips];

	Double_t ex, ey, exdiff, eydiff;
	Int_t multx = 0; Int_t multy = 0;

	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("Pspx1_xPrecal*",1);
	fChain->SetBranchStatus("Pspx1_yPrecal*",1);
	fChain->SetBranchStatus("Pspx2_xPrecal*",1);
	fChain->SetBranchStatus("Pspx2_yPrecal*",1);
	fChain->SetBranchStatus("Pspx3_xPrecal*",1);
	fChain->SetBranchStatus("Pspx3_yPrecal*",1);

	for (Long64_t jentry=0; jentry<neve*nentries;jentry++){
		Long64_t ientry = fChain->LoadTree(jentry);
		if (ientry < 0) break;

		nb = fChain->GetEntry(jentry);
		//cout << nb<<" "<< ientry << " " << jentry<< endl;
		nbytes += nb;
		if(nb==4) continue;
		if(nb==0) {cout<< "Something went terribly wrong!" << endl; break;};

		if (jentry%verbosity_step_size==0){
			cout << jentry << " of " << nentries << " == "<< 100. * jentry / nentries <<"% of the events analyzed\r"<< flush;
		}

		//	Reseting all values
		multx = Pspx1_xPrecal_;
		multy = Pspx1_yPrecal_;
		Int_t time_x = 0; Int_t time_y = 0;
		Int_t timediff_x = 0; Int_t timediff_y = 0;
		ex = 0; ey = 0;
		exdiff = 0; eydiff = 0;
		
		//timing calibratiom
        Int_t time_y0 = Pspx1_yPrecal_fTime1[0];
        Int_t time_y1 = Pspx1_yPrecal_fTime2[0];
        Int_t time_x0 = Pspx1_xPrecal_fTime1[0];
        Int_t time_x1 = Pspx1_xPrecal_fTime2[0]; 

	Int_t trigger_x = Pspx1_xPrecal_fTrigger[0];
	Int_t trigger_y = Pspx1_yPrecal_fTrigger[0];

        timediff_x = time_x0 - time_x1;
        time_x = time_x0 + time_x1 - 2.*trigger_x;
        timediff_y = time_y0 - time_y1;
        time_y = time_y0 + time_y1 - 2.*trigger_y;

        hTimeX[Pspx1_xPrecal_fStrip[0]-1]->Fill((time_x0+time_x1-2.*trigger_x)/2.);
        hTimeY[Pspx1_yPrecal_fStrip[0]-1]->Fill((time_y0+time_y1-2.*trigger_y)/2.);
        hTimeDiffX[Pspx1_xPrecal_fStrip[0]-1]->Fill(timediff_x/time_x);
        hTimeDiffY[Pspx1_yPrecal_fStrip[0]-1]->Fill(timediff_y/time_y);

		// cout << multx << " " << multy << endl;
		// middle interstrip event
		if (multx == 1 && multy == 2){
			Int_t str0 = Pspx1_yPrecal_fStrip[0];
			Int_t str1 = Pspx1_yPrecal_fStrip[1];
			if((str0 == strip_b0 && str1 == strip_b1) || (str1 == strip_b0 && str0 == strip_b1)){ 
				ex = Pspx1_xPrecal_fEnergy1[0] + Pspx1_xPrecal_fEnergy2[0]; 
				exdiff = Pspx1_xPrecal_fEnergy1[0] - Pspx1_xPrecal_fEnergy2[0]; 
				ey = Pspx1_yPrecal_fEnergy1[0] + Pspx1_yPrecal_fEnergy2[0]/* + Pspx1_yPrecal_fEnergy1[1] + Pspx1_yPrecal_fEnergy2[1]*/; 
				// if(((ex < ey + rangeE && ex > ey - rangeE) || (ey < ex + rangeE && ey > ex - rangeE))) {
					hPosX[Pspx1_xPrecal_fStrip[0] - 1]->Fill(exdiff/ex);
					counter_front++;
				// }
			}
		}

		// middle interstrip event
		if (multx == 2 && multy == 1){
			Int_t str0 = Pspx1_xPrecal_fStrip[0];
			Int_t str1 = Pspx1_xPrecal_fStrip[1];
			if((str0 == strip_f0 && str1 == strip_f1) || (str1 == strip_f0 && str0 == strip_f1)){ 
				ex = Pspx1_xPrecal_fEnergy1[0] + Pspx1_xPrecal_fEnergy2[0]/* + Pspx1_xPrecal_fEnergy1[1] + Pspx1_xPrecal_fEnergy2[1]*/; 
				ey = Pspx1_yPrecal_fEnergy1[0] + Pspx1_yPrecal_fEnergy2[0]; 
				eydiff = Pspx1_yPrecal_fEnergy1[0] - Pspx1_yPrecal_fEnergy2[0]; 
				// if(((ex < ey + rangeE && ex > ey - rangeE) || (ey < ex + rangeE && ey > ex - rangeE))) {
					hPosY[Pspx1_yPrecal_fStrip[0] - 1]->Fill(eydiff/ey);
					counter_back++;
				// }
			}
		}
		if((multx == 2 && multy == 1) || (multx == 1 && multy == 2)){
			hinterstrip->Fill(
				(Pspx1_xPrecal_fEnergy1[0] - Pspx1_xPrecal_fEnergy2[0])/(Pspx1_xPrecal_fEnergy1[0] + Pspx1_xPrecal_fEnergy2[0]),
				(Pspx1_yPrecal_fEnergy1[0] - Pspx1_yPrecal_fEnergy2[0])/(Pspx1_yPrecal_fEnergy1[0] + Pspx1_yPrecal_fEnergy2[0]));
		}
		hstrip->Fill(Pspx1_xPrecal_fStrip[0],Pspx1_yPrecal_fStrip[0]);
	}

	cout << "\nCounter front: " << counter_front << endl;
	cout << "Counter back: " << counter_back << endl;

	// --- FIT ---
	TF1* gaus_front = new TF1("gaus_front", "gaus(0)+pol0(3)", -0.06, 0.06);
	TF1* gaus_back = new TF1("gaus_back", "gaus(0)+pol0(3)", -0.06, 0.06);
	Double_t fitrange = 0.01;

	Double_t mean_front_middle[nStrips];
	Double_t sigma_front_middle[nStrips];

	Double_t mean_back_middle[nStrips];
	Double_t sigma_back_middle[nStrips];

	TSpectrum* s = new TSpectrum(2);
	Int_t nfound = 0;
	Double_t* xpeaks = nullptr;

	for (UInt_t i = 0; i <nStrips; i++) {
		cout << "Fit Strip " << i << endl;

		// Front
		nfound = s->Search(hPosX[i], 2, "", 0.20);
		xpeaks = s->GetPositionX();

		gaus_front->SetParameters(20, xpeaks[0], 0.02, 5);
		gaus_front->SetParLimits(0, 1, 5000);
		gaus_front->SetParLimits(1, xpeaks[0] - 0.01, xpeaks[0] + 0.01);
		gaus_front->SetParLimits(2, 0.0001, 0.2);
		// gaus_front->SetParLimits(3, 0, 600);
		gaus_front->FixParameter(3, 0);

		hPosX[i]->Fit("gaus_front", "IQ+", "", xpeaks[0] - fitrange, xpeaks[0] + fitrange);
		mean_front_middle[i] = gaus_front->GetParameter(1);

		// Back
		nfound = s->Search(hPosY[i], 2, "", 0.20);
		xpeaks = s->GetPositionX();

		gaus_back->SetParameters(20, xpeaks[0], 0.02, 5);
		gaus_back->SetParLimits(0, 1, 5000);
		gaus_back->SetParLimits(1, xpeaks[0] - 0.01, xpeaks[0] + 0.01);
		gaus_back->SetParLimits(2, 0.0001, 0.2);
		// gaus_back->SetParLimits(3, 0, 600);
		gaus_back->FixParameter(3, 0);

		hPosY[i]->Fit("gaus_back", "IQ+", "", xpeaks[0] - fitrange, xpeaks[0] + fitrange);

		mean_back_middle[i] = gaus_back->GetParameter(1);

		//if(i==17)
		//    gaus_front->SaveAs("/home/land/PSPX/plots/test.root");
		cout << i+1 << ", " << mean_front_middle[i] << ", " << mean_back_middle[i] << endl;
		//cout << i+1 << ", " << PspxPrecalData_fEnergy1[i] << ", " <<  PspxPrecalData_fEnergy2[i] << endl;
	}


	// --- Outfile ---
	if(save && !check){
		TFile * f_precal_histos = new TFile(Form("delme_pspx_precal_histos_bcal_run%i_to_run%i.root",run_min,run_max),"RECREATE");

		f_precal_histos->cd();
		for (UInt_t k = 0; k < nStrips; k++) {
			hPosX[k]->Write();
			hPosY[k]->Write();
		}
	} else if(save && check){
		TFile * f_precal_histos = new TFile(Form("delme_pspx_precal_histos_acal_run%i_to_run%i.root",run_min,run_max),"RECREATE");

		f_precal_histos->cd();
		for (UInt_t k = 0; k < nStrips; k++) { 
			hPosX[k]->Write();
			hPosY[k]->Write();
		}
	}


	// --- CALIB PARAM ---
	// Front
	cout << "1.0 1.0 32.0 \\" << endl;
	for (UInt_t i = 0; i < nStrips; i++) {
		cout << "  " << i + 1 << ".0 " << float(strip_b0)/float(32-strip_b0)*(1 + mean_front_middle[i]) / (1 - mean_front_middle[i]) << " 0.0 0.0 \\" << endl;
	}
	// Back
	cout << "1.0 2.0 32.0 \\" << endl;
	for (UInt_t i = 0; i < nStrips; i++) {
		cout << "  " << i + 1 << ".0 " << float(32-strip_f0)/float(strip_f0)*(1 + mean_back_middle[i]) / (1 - mean_back_middle[i]) << " 0.0 0.0 \\" << endl;
	}

	// --- PLOTS ---
	gStyle->SetOptStat(110011);
	gStyle->SetOptFit(0011);
	gStyle->SetPalette(1);

	for (UInt_t k = 0; k < nStrips; k++) {
		hPosX[k]->SetLineColor(4);
		hPosY[k]->SetLineColor(4);
	}
	TCanvas* cx1 = new TCanvas("PSPX1_X_center_part1", "PSPX1 X Center Part 1", 1400, 700);
	cx1->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		cx1->cd(1 + m);
		hPosX[8 + m]->Draw();
	}
	
	TCanvas* cx2 = new TCanvas("PSPX1_X_center_part2", "PSPX1 X Center Part 2", 1400, 700);
	cx2->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		cx2->cd(1 + m);
		hPosX[16 + m]->Draw();
	}
	
	TCanvas* cy1 = new TCanvas("PSPX1_Y_center_part1", "PSPX1 Y Center Part 1", 1400, 700);
	cy1->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		cy1->cd(1 + m);
		hPosY[8 + m]->Draw();
	}

	TCanvas* cy2 = new TCanvas("PSPX1_Y_center_part2", "PSPX1 Y Center Part 2", 1400, 700);
	cy2->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		cy2->cd(1 + m);
		hPosY[16 + m]->Draw();
	}

	TCanvas* ct1 = new TCanvas("PSPX1_X_time_part2", "PSPX1 X Time Part 2", 1400, 700);
	ct1->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		ct1->cd(1 + m);
		hTimeX[16 + m]->Draw();
	}

	TCanvas* ct2 = new TCanvas("PSPX1_Y_time_part2", "PSPX1 Y Time Part 2", 1400, 700);
	ct2->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		ct2->cd(1 + m);
		hTimeY[16 + m]->Draw();
	}

	TCanvas* ctd1 = new TCanvas("PSPX1_X_timediff_part2", "PSPX1 X TimeDiff Part 2", 1400, 700);
	ctd1->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		ctd1->cd(1 + m);
		hTimeDiffX[16 + m]->Draw();
	}

	TCanvas* ctd2 = new TCanvas("PSPX1_Y_timediff_part2", "PSPX1 Y TimeDiff Part 2", 1400, 700);
	ctd2->Divide(4, 2);
	for (Int_t m = 0; m < 8; m++) {
		ctd2->cd(1 + m);
		hTimeDiffY[16 + m]->Draw();
	}

	TCanvas* ccheck = new TCanvas("ccheck","Check",1400,700);
	ccheck->Divide(2,1);
	ccheck->cd(1);
	hinterstrip->Draw("colz");
	ccheck->cd(2);
	hstrip->Draw("colz");
}
