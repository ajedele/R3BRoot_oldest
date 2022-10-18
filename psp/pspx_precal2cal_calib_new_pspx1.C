/*
 * Automated calibration for PSPX1 detectors.
 *
 * HOWTO use: (not yet)
 * run script to get Mapped2Precal calibration parameters ("root -l pspx_mapped2precal_calib_frontback.C+(true)")
 * and copy them into the corresponding parameter file
 *
 * WARNING: The arrays in this script start counting with 0 (i.e. 0-4 for detectors, 0-15 for strips, 0-64 for
 * channels), but
 * the printed output is for an array gain starting with 1
 */

#include "TApplication.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TObject.h"
#include "TStyle.h"
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"

#include "TF1.h"
#include "TLatex.h"

#include <iostream>

#include "/u/sstorck/R3BRoot/r3bbase/R3BEventHeader.h"
#include "FairEventHeader.h"

using namespace std;
using namespace TMath;

void pspx_precal2cal_calib_new_pspx1(UInt_t run_min = 200,
                                     UInt_t run_max = 200,
                                     bool save = false,
                                     Double_t maxE = 10000,
                                     Double_t minE = 0,
                                     Double_t threshold = 10000,
                                     Double_t setE = 5000)
{

    UInt_t verbosity_step_size = 1000000; // number, how often progress of skript is reported

    // Detector specific variables
    UInt_t nStrips = 32;
    UInt_t npeaks = nStrips - 1;
    UInt_t run = 1;
    Double_t neve = 1.;

    // Input files
    TChain* fChain = new TChain("evt");
    TString filename;
    for (UInt_t i = run_min; i <= run_max; i++)
    {
        filename = Form("/lustre/land/sstorck/rootfiles/s515/main_run%03d_r3broot_precal.root", i);
        cout << "Open " << filename << endl;
        fChain->Add(filename);
    }

#include "new_class.h"
    fChain->SetMakeClass(1);

    // Calibration parameters
    Double_t mean_strip[nStrips];
    Double_t gain_strip[nStrips];

    for (UInt_t i = 0; i < nStrips; i++)
    {
        mean_strip[i] = 0;
        gain_strip[i] = 0;
    }

    // --- HISTOGRAMS ---
    Double_t minPos = -0.5;
    Double_t maxPos = 0.5;

    Double_t binPos = 200;

    Double_t binE = 500;

    Double_t titlesize = 0.045;

    TH2F* H_test_corr = new TH2F("H_test_corr", "H_test_corr", 500, minE, maxE, 500, minE, maxE);
    Double_t dummy1 = 0;
    Double_t dummy2 = 0;

    TH1F* hEx = new TH1F("energy_x", "Energy sum, X", binE, minE, maxE);
    TH1F* hEy = new TH1F("energy_y", "Energy sum, Y", binE, minE, maxE);
    TH1F* hMultx = new TH1F("mult_x", "Multiplicity, Front", 10, 0, 10);
    TH1F* hMulty = new TH1F("mult_y", "Multiplicity, Back", 10, 0, 10);
    TH1F* hStripsx = new TH1F("strips_x", "Strips, x Axis", nStrips, 1, nStrips + 1);
    TH1F* hStripsy = new TH1F("strips_y", "Strips, y Axis", nStrips, 1, nStrips + 1);

    TH1F* hEx_mult1 = new TH1F("energy_front_mult1", "Energy sum, X, Mult=1", binE, minE, maxE);
    TH1F* hEy_mult1 = new TH1F("energy_back_mult1", "Energy sum, Y, Mult=1", binE, minE, maxE);

    hEx_mult1->GetXaxis()->SetTitle("Energy/arb.u.");
    hEx_mult1->GetXaxis()->SetTitleSize(titlesize);
    hEx_mult1->GetYaxis()->SetTitle("Counts");
    hEx_mult1->GetYaxis()->SetTitleSize(titlesize);
    hEx_mult1->GetYaxis()->SetTitleOffset(1);

    hEy_mult1->GetXaxis()->SetTitle("Energy/arb.u.");
    hEy_mult1->GetXaxis()->SetTitleSize(titlesize);
    hEy_mult1->GetYaxis()->SetTitle("Counts");
    hEy_mult1->GetYaxis()->SetTitleSize(titlesize);
    hEy_mult1->GetYaxis()->SetTitleOffset(1);

    TH1F* hE_stripX[nStrips];
    TH1F* hE_stripY[nStrips];

    for (UInt_t k = 0; k < nStrips; k++)
    {
        hE_stripX[k] = new TH1F(Form("ene_strip_front%d", k), Form("Energy, Strip %d, Front", k + 1), binE, minE, maxE);
        hE_stripY[k] =
            new TH1F(Form("ene_strip_back%d", k), Form("Energy, Strip %d, Back", k + 1 + nStrips), binE, minE, maxE);

        hE_stripX[k]->GetXaxis()->SetTitle("Energy/arb.u.");
        hE_stripX[k]->GetXaxis()->SetTitleSize(titlesize);
        hE_stripX[k]->GetYaxis()->SetTitle("Counts");
        hE_stripX[k]->GetYaxis()->SetTitleSize(titlesize);
        hE_stripX[k]->GetYaxis()->SetTitleOffset(1);

        hE_stripY[k]->GetXaxis()->SetTitle("Energy/arb.u.");
        hE_stripY[k]->GetXaxis()->SetTitleSize(titlesize);
        hE_stripY[k]->GetYaxis()->SetTitle("Counts");
        hE_stripY[k]->GetYaxis()->SetTitleSize(titlesize);
        hE_stripY[k]->GetYaxis()->SetTitleOffset(1);
    }

    // --- EVENTLOOP ---
    if (fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    Int_t mult_front, mult_back;

    Int_t counter_front = 0;
    Int_t counter_back = 0;

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("Pspx1_xPrecal*", 1);
    fChain->SetBranchStatus("Pspx1_yPrecal*", 1);

    for (Long64_t jentry = 0; jentry < neve * nentries; jentry++)
    {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0)
            break;

        nb = fChain->GetEntry(jentry);
        // cout << nb<<" "<< ientry << " " << jentry<< endl;
        nbytes += nb;
        if (nb == 4)
            continue;
        if (nb == 0)
        {
            cout << "Something went terribly wrong!" << endl;
            break;
        };

        if (jentry % verbosity_step_size == 0)
        {
            cout << jentry << " of " << nentries << " == " << 100. * jentry / nentries << "% of the events analyzed\r"
                 << flush;
        }

        Int_t multx = Pspx1_xPrecal_;
        Int_t multy = Pspx1_yPrecal_;

        hMultx->Fill(multx);
        hMulty->Fill(multy);
        hEx->Fill(Pspx1_xPrecal_fEnergy1[0] + Pspx1_xPrecal_fEnergy2[0]);
        hEy->Fill(Pspx1_yPrecal_fEnergy1[0] + Pspx1_yPrecal_fEnergy2[0]);

        if (multx == 2 && multy == 1)
            counter_back++;
        else if (multx == 1 && multy == 2)
            counter_front++;

        // calculate total energy for each strip
        if (multx == 1)
        {
            for (Int_t i = 0; i < Pspx1_xPrecal_; i++)
            {
                if (!isnan(Pspx1_xPrecal_fEnergy1[i]) && !isnan(Pspx1_xPrecal_fEnergy2[i]))
                {
                    hE_stripX[Pspx1_xPrecal_fStrip[i] - 1]->Fill(Pspx1_xPrecal_fEnergy1[i] + Pspx1_xPrecal_fEnergy2[i]);
                    hEx_mult1->Fill(Pspx1_xPrecal_fEnergy1[i] + Pspx1_xPrecal_fEnergy2[i]);
                }
            }
        }
        if (multy == 1)
        {
            for (Int_t i = 0; i < Pspx1_yPrecal_; i++)
            {
                if (!isnan(Pspx1_yPrecal_fEnergy1[i]) && !isnan(Pspx1_yPrecal_fEnergy2[i]))
                {
                    hE_stripY[Pspx1_yPrecal_fStrip[i] - 1]->Fill(Pspx1_yPrecal_fEnergy1[i] + Pspx1_yPrecal_fEnergy2[i]);
                    hEy_mult1->Fill(Pspx1_yPrecal_fEnergy1[i] + Pspx1_yPrecal_fEnergy2[i]);
                }
            }
        }

        // Test correlation
        if (multx == 1 && multy == 1)
        {
            dummy1 = 0;
            dummy2 = 0;
            if (Pspx1_xPrecal_fStrip[0] == 17 && !isnan(Pspx1_xPrecal_fEnergy1[0]) && !isnan(Pspx1_xPrecal_fEnergy2[0]))
            {
                if (abs(Pspx1_xPrecal_fEnergy1[0]) > threshold && abs(Pspx1_xPrecal_fEnergy2[0]) > threshold)
                {
                    dummy1 = Pspx1_xPrecal_fEnergy1[0] + Pspx1_xPrecal_fEnergy2[0];
                }
            }
            if (Pspx1_yPrecal_fStrip[0] == 12 && !isnan(Pspx1_yPrecal_fEnergy1[0]) && !isnan(Pspx1_yPrecal_fEnergy2[0]))
            {
                if (abs(Pspx1_yPrecal_fEnergy1[0]) > threshold && abs(Pspx1_yPrecal_fEnergy2[0]) > threshold)
                {
                    dummy2 = Pspx1_yPrecal_fEnergy1[0] + Pspx1_yPrecal_fEnergy2[0];
                }
            }
            if (dummy1 != 0 && dummy2 != 0)
            {
                H_test_corr->Fill(dummy1, dummy2);
            }
        }
    }

    cout << "\nCounter front: " << counter_front << endl;
    cout << "Counter back: " << counter_back << endl;

    // Analysing energy spectra, calibration and resolution
    Double_t const_strip_front[nStrips];
    Double_t const_strip_back[nStrips];
    Double_t mean_strip_front[nStrips];
    Double_t mean_strip_back[nStrips];
    Double_t sigma_strip_front[nStrips];
    Double_t sigma_strip_back[nStrips];

    Double_t resolution_front = 0;
    Double_t resolution_back = 0;

    for (UInt_t i = 0; i < nStrips; i++)
    {
        const_strip_front[i] = 0;
        const_strip_front[i] = 0;
        mean_strip_front[i] = 0;
        mean_strip_back[i] = 0;
        sigma_strip_front[i] = 0;
        sigma_strip_back[i] = 0;
    }

    TF1* gaus_Ex = new TF1("gaus_Ex", "gaus(0)", minE, maxE); // minE + 0.45 * (maxE - minE)
    TF1* gaus_Ey = new TF1("gaus_Ey", "gaus(0)", minE, maxE); //-minE - 0.45 * (maxE - minE)

    // test fit with spectrum
    TSpectrum* s = new TSpectrum(15);
    Int_t nfound = 0;
    Double_t* xpeaks;
    Double_t fit_range = 10000;
    Double_t fit_min = 0;
    Double_t fit_max = 0;

    for (UInt_t k = 0; k < nStrips; k++)
    {

        // test fit with spectrum
        nfound = s->Search(hE_stripX[k], 2, "", 0.60);
        xpeaks = s->GetPositionX();
        gaus_Ex->SetParameter(0, 5000);
        gaus_Ex->SetParameter(1, xpeaks[0]);
        gaus_Ex->SetParameter(2, 0.02 * xpeaks[0]);
        fit_min = xpeaks[0] - fit_range;
        fit_max = xpeaks[0] + fit_range;

        hE_stripX[k]->Fit("gaus_Ex", "BIQ", "", fit_min, fit_max);
        mean_strip_front[k] = gaus_Ex->GetParameter(1);
        sigma_strip_front[k] = gaus_Ex->GetParameter(2);
        resolution_front += abs(sigma_strip_front[k] / mean_strip_front[k]) * 100.;

        nfound = s->Search(hE_stripY[k], 2, "", 0.60);
        xpeaks = s->GetPositionX();
        gaus_Ey->SetParameter(0, 5000);
        gaus_Ey->SetParameter(1, xpeaks[0]);
        gaus_Ey->SetParameter(2, 0.02 * xpeaks[0]);
        fit_min = xpeaks[0] - fit_range;
        fit_max = xpeaks[0] + fit_range;

        hE_stripY[k]->Fit("gaus_Ey", "BIQ", "", fit_min, fit_max);
        mean_strip_back[k] = gaus_Ey->GetParameter(1);
        sigma_strip_back[k] = gaus_Ey->GetParameter(2);
        resolution_back += abs(sigma_strip_back[k] / mean_strip_back[k]) * 100.;
        // end of spectrum

        cout << "Strip " << k + 1 << " Front: " << abs(sigma_strip_front[k] / mean_strip_front[k]) * 100.
             << " % Back: " << abs(sigma_strip_back[k] / mean_strip_back[k]) * 100. << " %" << endl;

        cout << k + 1 << ", front: " << mean_strip_front[k] << ", back: " << mean_strip_back[k] << endl;
    }
    cout << "Resolution Mean Front: " << resolution_front / nStrips << " % Back: " << resolution_back / nStrips << " %"
         << endl;

    gaus_Ex->SetParameter(0, 5000);
    gaus_Ex->SetParameter(1, maxE / 4.);
    gaus_Ex->SetParameter(2, maxE / 20.);

    gaus_Ey->SetParameter(0, 5000);
    gaus_Ey->SetParameter(1, maxE / 4.);
    gaus_Ey->SetParameter(2, maxE / 20.);

    hEx_mult1->Fit("gaus_Ex", "IQR"); // 0 to make fit invisible
    Double_t const_front_mult1 = gaus_Ex->GetParameter(0);
    Double_t mean_front_mult1 = gaus_Ex->GetParameter(1);
    Double_t sigma_front_mult1 = gaus_Ex->GetParameter(2);

    hEy_mult1->Fit("gaus_Ey", "IQR");
    Double_t const_back_mult1 = gaus_Ey->GetParameter(0);
    Double_t mean_back_mult1 = gaus_Ey->GetParameter(1);
    Double_t sigma_back_mult1 = gaus_Ey->GetParameter(2);

    cout << "Resolution Front: " << abs(sigma_front_mult1 / mean_front_mult1) * 100.
         << " % Back: " << abs(sigma_back_mult1 / mean_back_mult1) * 100. << " %" << endl;

    // --- CALIB PARAM ---

    cout << "Calib Parameter " << endl;
    // Front
    cout << "1.0 1.0 32.0 \\" << endl;
    for (UInt_t i = 0; i < nStrips; i++)
    {
        cout << "  " << i + 1 << ".0 " << setE / mean_strip_front[i] << " \\" << endl;
    }
    // Back
    cout << "1.0 2.0 32.0 \\" << endl;
    for (UInt_t i = 0; i < nStrips; i++)
    {
        cout << "  " << i + 1 << ".0 " << setE / mean_strip_back[i] << " \\" << endl;
    }

    // --- Outfile ---
    if (save)
    {
        TFile* f_cal_histos =
            new TFile(Form("../plots/delme_pspx_cal_histos_bcal_%i_det1_run%i_to_run%i.root", 181126, run_min, run_max),
                      "RECREATE");

        f_cal_histos->cd();
        for (UInt_t k = 0; k < nStrips; k++)
        {
            hE_stripX[k]->Write();
            hE_stripY[k]->Write();
        }

        hEx->Write();
        hEy->Write();
        hMultx->Write();
        hMulty->Write();
        hStripsx->Write();
        hStripsy->Write();
        hEx_mult1->Write();
        hEy_mult1->Write();
    }

    // --- PLOTS ---
    gStyle->SetOptStat(111111);
    // gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas* c1 = new TCanvas("PSPX1_x_center_ene_part1", "PSPX1 X Center Energy Part1", 1400, 700);
    c1->Divide(4, 4);

    for (Int_t m = 0; m < 16; m++)
    {
        c1->cd(1 + m);
        hE_stripX[m]->Draw();
    }

    TCanvas* c2 = new TCanvas("PSPX1_x_center_ene_part2", "PSPX1 X Center Energy Part2", 1400, 700);
    c2->Divide(4, 4);

    for (Int_t m = 0; m < 16; m++)
    {
        c2->cd(1 + m);
        hE_stripX[16 + m]->Draw();
    }

    TCanvas* c4 = new TCanvas("PSPX1_y_center_ene_part1", "PSPX1 Y Center Energy Part 1", 1400, 700);
    c4->Divide(4, 4);

    for (Int_t m = 0; m < 16; m++)
    {
        c4->cd(1 + m);
        hE_stripY[m]->Draw();
    }

    TCanvas* c5 = new TCanvas("PSPX1_y_center_ene_part2", "PSPX1 Y Center Energy Part 2", 1400, 700);
    c5->Divide(4, 4);

    for (Int_t m = 0; m < 16; m++)
    {
        c5->cd(1 + m);
        hE_stripY[16 + m]->Draw();
    }

    TCanvas* c6 = new TCanvas("PSPX1_enery", "PSPX1 Energy", 1400, 700);
    c6->Divide(2, 1);

    c6->cd(1);
    hEx_mult1->Draw();

    TLatex res_front;
    res_front.DrawLatex(minE + 0.1 * (maxE - minE),
                        0.5 * hEx_mult1->GetMaximum(),
                        Form("#sigma_{front} = %.2f %%", abs(sigma_front_mult1 / mean_front_mult1) * 100.));

    c6->cd(2);
    hEy_mult1->Draw();

    TLatex res_back;
    res_back.DrawLatex(minE + 0.1 * (maxE - minE),
                       0.5 * hEy_mult1->GetMaximum(),
                       Form("#sigma_{back} = %.2f %%", abs(sigma_back_mult1 / mean_back_mult1) * 100.));

    Double_t aspect_ratio = TMath::Sqrt(2);

    TCanvas* c34 = new TCanvas("c34", "Test Corr", 1000, 1000);
    c34->cd();
    H_test_corr->Draw("colz");
}
