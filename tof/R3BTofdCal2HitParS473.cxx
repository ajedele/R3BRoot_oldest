/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ----------------------------------------------------------------
// -----            R3BTofdCal2HitPar             -----
// -----           Created May 2016 by M.Heil             -----
// ----------------------------------------------------------------

/* Some notes:
 *
 *
 */

#include "R3BTofdCal2HitParS473.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include "TLine.h"
#include "TLegend.h"

#include <iostream>
#include <stdlib.h>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdCal2HitParS473::R3BTofdCal2HitParS473()
    : FairTask("R3BTofdCal2HitParS473", 1)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofdY(0.)
    , fTofdQ(0.)
    , fParaFile("")
    , maxEvents(0)
    , fTofdToTLow(10.)
    , fTofdToTHigh(150.)
    , fLastEvent(0)
    , fMultCut(false)
    , fTofdSmiley(false)
    , fTofdWalk(false)
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            fhTotPm1[i][j] = NULL;
            fhTotPm2[i][j] = NULL;
            fhTot1vsTot2[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
            fhSqrtQvsPosToT[i][j] = NULL;
        }
    }
}

R3BTofdCal2HitParS473::R3BTofdCal2HitParS473(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofdY(0.)
    , fTofdQ(0.)
    , fParaFile("")
    , maxEvents(0)
    , fTofdToTLow(10.)
    , fTofdToTHigh(150.)
    , fLastEvent(0)
    , fMultCut(false)
    , fTofdSmiley(false)
    , fTofdWalk(false)
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            fhTotPm1[i][j] = NULL;
            fhTotPm2[i][j] = NULL;
            fhTot1vsTot2[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
	    fhSqrtQvsPosToT[i][j] = NULL;
        }
    }
}

R3BTofdCal2HitParS473::~R3BTofdCal2HitParS473()
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        if (fhTdiff[i])
            delete fhTdiff[i];
        if (fhTsync[i])
            delete fhTsync[i];
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fhTotPm1[i][j])
                delete fhTotPm1[i][j];
            if (fhTotPm2[i][j])
                delete fhTotPm2[i][j];
            if (fhTot1vsTot2[i][j])
                delete fhTot1vsTot2[i][j];
            if (fhTot1vsPos[i][j])
                delete fhTot1vsPos[i][j];
            if (fhTot2vsPos[i][j])
                delete fhTot2vsPos[i][j];
	    if(fhSqrtQvsPosToT[i][j])
		delete fhSqrtQvsPosToT[i][j];
        }
    }
    if (fCal_Par)
    {
        delete fCal_Par;
    }
    if (fhToTvsDist){delete fhToTvsDist;}
}

InitStatus R3BTofdCal2HitParS473::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = (R3BEventHeader*)rm->GetObject("R3BEventHeader");
    // may be = NULL!

    fCalData = (TClonesArray*)rm->GetObject("TofdCal");
    if (!fCalData)
    {
        return kFATAL;
    }

    if (!fNofModules)
    {
        LOG(ERROR) << "R3BTofdCal2HitPar::Init() Number of modules not set. ";
        return kFATAL;
    }
    fCalItemsLos = (TClonesArray*)rm->GetObject("LosCal");
    if (NULL == fCalItemsLos)
        LOG(fatal) << "Branch LosCal not found";

    maxEvents = rm->CheckMaxEventNo();

    return kSUCCESS;
}

void R3BTofdCal2HitParS473::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCal_Par = (R3BTofdHitPar*)FairRuntimeDb::instance()->getContainer("TofdHitPar");
    if (!fCal_Par)
    {
        LOG(ERROR) << "R3BTofdCal2HitPar::Init() Couldn't get handle on TofdHitPar. ";
    }
    //	    fCal_Par->setChanged();
}

void R3BTofdCal2HitParS473::Exec(Option_t* option)
{
    if (fNEvents / 10000. == (int) fNEvents / 10000){
        cout << "Events: " << fNEvents << " / " << maxEvents << " (" << (Int_t) (fNEvents * 100. / maxEvents) << "%)"
	<< "         \r" << flush;
    }
    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Double_t timeLos = 0;
    Double_t time_r_V = nan(""), time_t_V = nan(""), time_l_V = nan(""), time_b_V = nan(""), time_rt_V = nan(""),
             time_lt_V = nan(""), time_lb_V = nan(""), time_rb_V = nan("");

    // Los detector
    if (fCalItemsLos)
    {
        Int_t nHits = fCalItemsLos->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosCalData* calData = (R3BLosCalData*)fCalItemsLos->At(ihit);

            Int_t iDet = calData->GetDetector();
            // Int_t iCha=calData->GetChannel();

            if (!(IS_NAN(calData->GetTimeV_ns(5))))
                time_r_V = calData->GetTimeV_ns(5);
            if (!(IS_NAN(calData->GetTimeV_ns(7))))
                time_t_V = calData->GetTimeV_ns(7);
            if (!(IS_NAN(calData->GetTimeV_ns(1))))
                time_l_V = calData->GetTimeV_ns(1);
            if (!(IS_NAN(calData->GetTimeV_ns(3))))
                time_b_V = calData->GetTimeV_ns(3);
            if (!(IS_NAN(calData->GetTimeV_ns(6))))
                time_rt_V = calData->GetTimeV_ns(6);
            if (!(IS_NAN(calData->GetTimeV_ns(0))))
                time_lt_V = calData->GetTimeV_ns(0);
            if (!(IS_NAN(calData->GetTimeV_ns(2))))
                time_lb_V = calData->GetTimeV_ns(2);
            if (!(IS_NAN(calData->GetTimeV_ns(4))))
                time_rb_V = calData->GetTimeV_ns(4);
        }
        timeLos = (time_r_V + time_t_V + time_l_V + time_b_V + time_rt_V + time_lt_V + time_lb_V + time_rb_V) / 8.;
    }

    // ToFD detector

    Int_t nHits = fCalData->GetEntries();

    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = (R3BTofdCalData*)fCalData->At(ihit);
        size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

        // std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << ' '
        //          << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << '\n';

        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit);
    }

    struct Hit
    {
	R3BTofdCalData* top;
	R3BTofdCalData* bot;
	Double_t timeLos;

    };
    vector<Hit> hitVec;
    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
        // for (auto it2 = it->second.top.begin(); it->second.top.end() != it2; ++it2) {
        // std::cout << "Top: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
        // for (auto it2 = it->second.bot.begin(); it->second.bot.end() != it2; ++it2) {
        // std::cout << "Bot: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
        auto const& top_vec = it->second.top;
        auto const& bot_vec = it->second.bot;
        size_t top_i = 0;
        size_t bot_i = 0;
        for (; top_i < top_vec.size() && bot_i < bot_vec.size();)
        {
            auto top = top_vec.at(top_i);
            auto bot = bot_vec.at(bot_i);
            auto top_ns = top->GetTimeLeading_ns();
            auto bot_ns = bot->GetTimeLeading_ns();

            auto dt = top_ns - bot_ns;
            // Handle wrap-around.
            auto dt_mod = fmod(dt + c_range_ns, c_range_ns);
            if (dt < 0)
            {
                // We're only interested in the short time-differences, so we
                // want to move the upper part of the coarse counter range close
                // to the lower range, i.e. we cut the middle of the range and
                // glue zero and the largest values together.
                dt_mod -= c_range_ns;
            }
            // std::cout << top_i << ' ' << bot_i << ": " << top_ns << ' ' << bot_ns << " = " << dt << ' ' <<
            // std::abs(dt_mod) << '\n';
            if (std::abs(dt_mod) < c_bar_coincidence_ns)
            {
                // Hit!

                // std::cout << "Hit!\n";
		Hit hit;
		hit.top = top;
		hit.bot = bot;
		hit.timeLos = timeLos;
		hitVec.push_back(hit);


                ++top_i;
                ++bot_i;

                // Increment events
                fNEvents += 1;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
            }
            else
            {
                ++bot_i;
            }
        }
    }
    fLastEvent++;

    bool skipEvt = false;
    if(fMultCut){
	    for(Int_t i = 0; i < hitVec.size(); i++){
		    for(Int_t j = 0; j < i; j++){
		    	if(hitVec.at(i).top->GetDetectorId() == hitVec.at(j).top->GetDetectorId()) skipEvt = true; 
		    }
	    }
    }

    if(!skipEvt){
	  for(Int_t i = 0; i < hitVec.size(); i++){

		  Hit hit = hitVec.at(i);
		  auto top = hit.top;
		  auto bot = hit.bot;

		  Int_t iPlane = top->GetDetectorId(); // 1..n
		  Int_t iBar = top->GetBarId();        // 1..n
		  if (iPlane > fNofPlanes)             // this also errors for iDetector==0
		  {
			  LOG(ERROR) << "R3BTofdCal2HitPar::Exec() : more detectors than expected! Det: " << iPlane
				  << " allowed are 1.." << fNofPlanes;
			  continue;
		  }
		  if (iBar > fPaddlesPerPlane) // same here
		  {
			  LOG(ERROR) << "R3BTofdCal2HitPar::Exec() : more bars then expected! Det: " << iBar
				  << " allowed are 1.." << fPaddlesPerPlane;
			  continue;
		  }

		  auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
		  auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);


		  auto top_ns = top->GetTimeLeading_ns();
		  auto bot_ns = bot->GetTimeLeading_ns();

		  // walk corrections
		  /*
		  bot_ns = bot_ns - walk(bot_tot);
		  top_ns = top_ns - walk(top_tot);
		  */

		  auto tdiff = bot_ns - top_ns;

		  // std::cout << "Saved: " << top->GetDetectorId() << ' ' << top->GetBarId() << ' '
		  //<< bot->GetDetectorId() << ' ' << bot->GetBarId() << ' '
		  //<< (top_ns + bot_ns)/2 << ' ' << tdiff << ' ' << top_tot << ' ' << bot_tot << '\n';

		  // create histograms
		  CreateHistograms(iPlane, iBar);

		  // fill histograms

		  fhTotPm1[iPlane - 1][iBar - 1]->Fill(bot_tot);
		  fhTotPm2[iPlane - 1][iBar - 1]->Fill(top_tot);
		  fhTot1vsTot2[iPlane - 1][iBar - 1]->Fill(top_tot, bot_tot);

		  // Time differences of one paddle **************************
		  fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

		  //  ToT vs Energy (== SQRT(tot1 * tot2)
		  auto posToT = 0.;
		  if (fTofdY == 0.)
			  posToT = TMath::Log(top_tot / bot_tot);
		  else
		  {
			  R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
			  if (!par)
			  {
				  LOG(ERROR) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
					  << ", Bar: " << top->GetBarId();
				  continue;
			  }
			  Double_t lambda = 1.;
			  if(fTofdQ != 0){lambda = par->GetLambda();}
			  posToT = lambda * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));
		  }
		  fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, sqrt(top_tot * bot_tot));


		  // ToF
		  auto ToF = (top_ns + bot_ns) / 2 - timeLos;
		  while (ToF < -c_range_ns / 2)
			  ToF += c_range_ns;
		  while (ToF > c_range_ns / 2)
			  ToF -= c_range_ns;

		  fhTsync[iPlane - 1]->Fill(iBar, ToF);

		  // cout << "test " << (top_ns + bot_ns)/2 <<"  " << timeLos << "  " << ToF << endl;

		  if (fTofdQ != 0)
		  {
			  R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
			  if (!par)
			  {
				  LOG(INFO) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
					  << ", Bar: " << top->GetBarId();
				  continue;
			  }
			  Double_t Pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

			  fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(Pos, bot_tot);
			  fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(Pos, top_tot);
			  
			  Double_t dist = 50.- Pos;
			  fhToTvsDist->Fill(dist, top_tot);
			  dist = -(-50. - Pos);
			  fhToTvsDist->Fill(dist, bot_tot);
		  }


	  }
   }
}


void R3BTofdCal2HitParS473::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 60.;

    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 4000, -20., 20.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }
    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, "", 50, 0, 50, 10000, -10, 90.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("ToF in ns");
    }

    if (NULL == fhTotPm1[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_Plane_%d_Bar_%d_PM_1", iPlane, iBar);
        fhTotPm1[iPlane - 1][iBar - 1] = new TH1F(strName, "", 300, 0., 300.);
        fhTotPm1[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM1 in ns");
        fhTotPm1[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Counts");
    }
    if (NULL == fhTotPm2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_Plane_%d_Bar_%d_PM_2", iPlane, iBar);
        fhTotPm2[iPlane - 1][iBar - 1] = new TH1F(strName, "", 300, 0., 300.);
        fhTotPm2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTotPm2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Counts");
    }
    if (NULL == fhTot1vsTot2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Plane_%d_Bar_%d_ToT1vsToT2", iPlane, iBar);
        fhTot1vsTot2[iPlane - 1][iBar - 1] = new TH2F(strName, "", 300, 0., 300., 300, 0., 300.);
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Tot of PM1 in ns");
    }
    if (NULL == fhTot1vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot1_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 400, 0., 200.);
        if (iPlane > 2)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 400, 0., 200.);
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Tot of PM1 in ns");
    }
    if (NULL == fhTot2vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot2_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 150, 50., 200.);
        if (iPlane > 2)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 150, 50., 200.);
    }
    if (NULL == fhSqrtQvsPosToT[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "SqrtQ_vs_PosToT_Plane_%d_Bar_%d", iPlane, iBar);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1] =
            new TH2F(strName, "", 20000, -100, 100, max_charge * 4, 0., max_charge * 4);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("sqrt(PM1*PM2)");
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position from ToT in cm");
    }
    if (NULL == fhToTvsDist)
    {
    	char strName[255];
	sprintf(strName, "ToT_vs_Distance");
	fhToTvsDist = new TH2F(strName, "", 20000, 0., 100., max_charge * 4, 0., max_charge * 4);
	fhToTvsDist->GetXaxis()->SetTitle("Distance to PMT in cm");
	fhToTvsDist->GetYaxis()->SetTitle("ToT of PMT in ns");
    }

}

void R3BTofdCal2HitParS473::FinishEvent()
{
    if (fCalItemsLos)
    {
        fCalItemsLos->Clear();
    }
}

void R3BTofdCal2HitParS473::FinishTask()
{

    cout << fLastEvent << " / " << maxEvents << " analysed for this calibration." << endl;

    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        if (fhTsync[i])
            fhTsync[i]->Write();
        if (fhTdiff[i])
            fhTdiff[i]->Write();
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fhTot1vsPos[i][j])
                fhTot1vsPos[i][j]->Write();
            if (fhTot2vsPos[i][j])
                fhTot2vsPos[i][j]->Write();
            if (fhTotPm1[i][j])
                fhTotPm1[i][j]->Write();
            if (fhTotPm2[i][j])
                fhTotPm2[i][j]->Write();
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write();
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write();
	    if(fhSqrtQvsPosToT[i][j])
	        fhSqrtQvsPosToT[i][j]->Write();
        }
    }

    // Determine time offset of the 2 PMTs of one paddle. This procedure
    // assumes a sweep run in the middle of the ToF wall horizontally.
    // Since all paddles are mounted vertically one can determine the offset.
    // Half of the offset is added to PM1 and half to PM2.

    if (fTofdY == 0 && fTofdQ == 0)
    {
        calcOffset();
        calcSync();
	calcToTOffset(fTofdToTLow, fTofdToTHigh);
    }
    else if (fTofdY != 0 && fTofdQ == 0)
    {
        calcVeff();
	calcLambda(fTofdToTLow, fTofdToTHigh);
    }

    Double_t para[4];
    Double_t min = -40;
    Double_t max = 40;
    if ((fTofdQ != 0) && !fTofdWalk)
    {
	setWalkPar();
	if(!fTofdSmiley){
            for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
            {
                for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
                {
                    if (fhTot1vsPos[i][j])
                    {
                        R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);

                        doubleExp(fhTot1vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar1a(para[0]);
                        par->SetPar1b(para[1]);
                        par->SetPar1c(para[2]);
                        par->SetPar1d(para[3]);
                    }

                    if (fhTot2vsPos[i][j])
                    {
                        R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                        doubleExp(fhTot2vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar2a(para[0]);
                        par->SetPar2b(para[1]);
                        par->SetPar2c(para[2]);
                        par->SetPar2d(para[3]);
                    }
                }
            }
        
           fCal_Par->setChanged();
	}else {
        
            LOG(WARNING) << "Calling function smiley";
            Double_t para2[4];
            Double_t min2 = -40.; // -40 effective bar length
            Double_t max2 = 40.;  // 40 effective bar length = 80 cm
                                  // we will use 50 here for some fit safety margin
            for (Int_t i = 0; i < fNofPlanes; i++)
            {
                for (Int_t j = 0; j < fPaddlesPerPlane; j++)
                {
                    if (fhSqrtQvsPosToT[i][j])
                    {
                        LOG(WARNING) << "Calling Plane " << i + 1 << " Bar " << j + 1;
                        R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                        smiley(fhSqrtQvsPosToT[i][j],
                               min2,
                               max2,
                               para2);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar1a(para2[0]);
                        par->SetPar1b(para2[1]);
                        par->SetPar1c(para2[2]);
                        par->SetPar1d(para2[3]);
                    }
                }
            }
            fCal_Par->setChanged();
	}
    }
    if((fTofdQ != 0) && fTofdWalk){
	Double_t para3[6];
    	fitWalk(fhToTvsDist,-40.,40.,para3);
	// TODO: set walk pars to fitted pars?
    }
}

void R3BTofdCal2HitParS473::calcOffset()
{
    Double_t sig[4][44] = { 0. };
    TCanvas* cOffset = new TCanvas("cOffset", "cOffset", 10, 10, 1000, 900);
    cOffset->Divide(2, 2);
    R3BTofdHitModulePar* mpar;
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        if (fhTdiff[i])
        {
            auto* h = (TH2F*) fhTdiff[i]->Clone();
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                mpar = new R3BTofdHitModulePar();
                Double_t offset = 0.;
                cOffset->cd(i + 1);
                h->Draw("colz");
                TH1F* histo_py = (TH1F*)h->ProjectionY("histo_py", j + 2, j + 2, "");
                histo_py->Rebin(4);
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.3, Max + 0.3);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1); // histo_py->GetXaxis()->GetBinCenter(binmax);
                LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " Offset  " << offset;
                mpar->SetPlane(i + 1);
                mpar->SetPaddle(j + 1);
                mpar->SetOffset1(-offset / 2.);
                mpar->SetOffset2(offset / 2.);
                fCal_Par->AddModulePar(mpar);
       	    }    
        }


    }

    fCal_Par->setChanged();
}

void R3BTofdCal2HitParS473::calcSync()
{
    TCanvas* cSync = new TCanvas("cSync", "cSync", 10, 10, 1000, 900);
    cSync->Divide(2, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        if (fhTsync[i])
        {
            auto* h = (TH2F*) fhTsync[i]->Clone();
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                cSync->cd(i + 1);
                h->Draw("colz");
                auto* histo_py = (TH1F*)h->ProjectionY("histo_py", j + 2, j + 2, "");
                R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                Double_t MaxEntry = histo_py->GetBinContent(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 10., Max + 10.);
                fgaus->SetParameters(MaxEntry, Max, 20);
                histo_py->Fit("fgaus", "QR0");
                Double_t sync = fgaus->GetParameter(1); // histo_py->GetXaxis()->GetBinCenter(binmax);
                par->SetSync(sync);
                LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " Sync  " << sync;
            }
        }
    }
    fCal_Par->setChanged();
}

void R3BTofdCal2HitParS473::calcToTOffset(Double_t totLow, Double_t totHigh)
{
	
    TCanvas* cToTOffset = new TCanvas("cToTOffset", "cToTOffset", 10, 10, 1000, 900);
    cToTOffset->Divide(1, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t offset = 0.;
            R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
	    if(fhSqrtQvsPosToT[i][j])
            {
                LOG(WARNING) << "Found histo SqrtQ_vs_PosToT_Plane_" << i + 1 << "_Bar_" << j + 1;
                TH2F* h = (TH2F*) fhSqrtQvsPosToT[i][j]->Clone();
                cToTOffset->cd(1);
                h->Draw("colz");
                auto* histo_py = (TH2F*)h->ProjectionX("histo_py", totLow, totHigh, "");
                cToTOffset->cd(2);
                histo_py->Rebin(2);
                histo_py->Draw();
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1(
                    "fgaus", "gaus(0)", Max - 0.2, Max + 0.2); // new TF1("fgaus", "gaus(0)", Max - 0.06, Max + 0.06);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1);
                fgaus->Draw("SAME");
                histo_py->SetAxisRange(Max - .5, Max + .5, "X");
                h->SetAxisRange(Max - .5, Max + .5, "X");
                h->SetAxisRange(totLow, totHigh, "Y");
                cToTOffset->Update();
                delete fgaus;
                delete h;
                delete histo_py;
            }
            LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " ToT Offset  " << offset << "\n";
            par->SetToTOffset1(sqrt(exp(offset)));
            par->SetToTOffset2(1. / sqrt(exp(offset)));
        }
    }
    fCal_Par->setChanged();
}

void R3BTofdCal2HitParS473::calcVeff()
{
    TCanvas* cVeff = new TCanvas("cVeff", "cVeff", 10, 10, 1000, 900);
    cVeff->Divide(2, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t max = 0.;
            Double_t veff = 7.;
            if (fhTdiff[i])
            {
                LOG(WARNING) << "Found histo Time_Diff_Plane_" << i + 1;
                R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                if (!par)
                {
                    LOG(INFO) << "Hit par not found, Plane: " << i + 1 << ", Bar: " << j + 1;
                    continue;
                }
                auto* h = (TH2F*) fhTdiff[i]->Clone();
                cVeff->cd(i + 1);
                h->Draw("colz");
                TH1F* histo_py = (TH1F*)h->ProjectionY("histo_py", j + 2, j + 2, "");
                Int_t binmax = histo_py->GetMaximumBin();
                max = histo_py->GetXaxis()->GetBinCenter(binmax);
                Double_t maxEntry = histo_py->GetBinContent(binmax);
                auto* fgaus = new TF1("fgaus", "gaus(0)", max - 0.3, max + 0.3); /// TODO: find best range
                fgaus->SetParameters(maxEntry, max, 20);
                histo_py->Fit("fgaus", "QR0");
                Double_t offset1 = par->GetOffset1();
                Double_t offset2 = par->GetOffset2();
                Double_t sync = par->GetSync();
                max = fgaus->GetParameter(1) + offset1 - offset2; /// TODO: needs to be tested
                // max = max+offset1-offset2;
                veff = fTofdY / max; // effective speed of light in [cm/s]
                LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " offset  " << par->GetOffset1();
                LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " max  " << max;
                LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " veff  " << veff;
                par->SetVeff(veff);
            }
        }
    }
    fCal_Par->setChanged();
}


void R3BTofdCal2HitParS473::calcLambda(Double_t totLow, Double_t totHigh)
{
    TCanvas* cToTOffset = new TCanvas("cLambda", "cLambda", 10, 10, 1000, 900);
    cToTOffset->Divide(1, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t offset = 0.;
            R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
            if (fhSqrtQvsPosToT[i][j])
            {
                LOG(WARNING) << "Found histo SqrtQ_vs_PosToT_Plane_" << i + 1 << "_Bar_" << j + 1;
                auto* h = (TH2F*) fhSqrtQvsPosToT[i][j]->Clone();
                cToTOffset->cd(1);
                h->Draw("colz");
                auto* histo_py = (TH2F*)h->ProjectionX("histo_py", totLow, totHigh, "");
                cToTOffset->cd(2);
                histo_py->Draw();
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.06, Max + 0.06);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1);
                fgaus->Draw("SAME");
                histo_py->SetAxisRange(Max - .5, Max + .5, "X");
                h->SetAxisRange(Max - .5, Max + .5, "X");
                h->SetAxisRange(totLow, totHigh, "Y");
                cToTOffset->Update();
                delete fgaus;
                delete h;
                delete histo_py;
            }
            else
                LOG(ERROR) << "Missing histo plane " << i + 1 << " bar " << j + 1;
            Double_t lambda = fTofdY / offset;
            LOG(WARNING) << " Plane  " << i + 1 << " Bar " << j + 1 << " ToT Offset  " << offset << " Lambda " << lambda
                         << "\n";
            par->SetLambda(lambda);
        }
    }
    fCal_Par->setChanged();
}

void R3BTofdCal2HitParS473::doubleExp(TH2F* histo, Double_t min, Double_t max, Double_t* para)
{
    // This fits the exponential decay of the light in a paddle. The 2 PMTs are fit with the same function but one
    // side will deliver negative attenuation parameters and the other positive.
    Double_t y[1000], x[1000];
    Int_t n = 0;
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = 0;
    }
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    TCanvas* cfit_exp = new TCanvas("cfit_exp", "fit exponential", 100, 100, 800, 800);
    cfit_exp->Clear();
    cfit_exp->Divide(1, 3);
    cfit_exp->cd(1);
    TH2F* histo1 = (TH2F*)histo->Clone();
    TH2F* histo2 = (TH2F*)histo->Clone();
    histo1->Draw("colz");
    cfit_exp->cd(2);
    for (Int_t i = 1; i < histo1->GetNbinsX() - 1; i++)
    {
        TH1F* histo_py = (TH1F*)histo1->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo1->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if ((x[n] < -40. || x[n] > 40.) || y[n] < 50.)
        {
            delete histo_py;
            continue;
        }
        if (histo_py->GetMaximum() > 5)
            n++;
        delete histo_py;
    }
    gr1 = new TGraph(n, x, y);
    gr1->Draw("A*");
    TF1* f1 = new TF1("f1", "[0]*(exp(-[1]*(x+100.))+exp(-[2]*(x+100.)))+[3]", min, max);
    f1->SetParameters(520., 0.001, 17234, -485.);
    f1->SetLineColor(2);
    gr1->Fit("f1", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f1->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    // fit again but with more information and better cuts
    n = 0;
    cfit_exp->cd(3);
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        Double_t pos = histo2->GetXaxis()->GetBinCenter(i);
        Double_t ymean = para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3];
        histo2->SetAxisRange(ymean - 5., ymean + 5., "Y");
        histo2->Draw("colz");
        TH1F* histo_py = (TH1F*)histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 2)
            n++;
        delete histo_py;
    }
    gr2 = new TGraph(n, x, y);
    gr2->Draw("A*");
    TF1* f2 = new TF1("f2", "[0]*(exp(-[1]*(x+100.))+exp(-[2]*(x+100.)))+[3]", min, max);
    f2->SetParameters(para[0], para[1], para[2], para[3]);
    f2->SetLineColor(2);
    gr2->Fit("f2", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f2->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    cfit_exp->Update();
    // gPad->WaitPrimitive();
    gSystem->Sleep(3000);
    delete gr1;
    delete gr2;
    delete f1;
    delete f2;
}

void R3BTofdCal2HitParS473::smiley(TH2F* histo, Double_t min, Double_t max, Double_t* para)
{
    // This fits the smiley: Sqrt(q1*q2) returns position dependent charge, we fit that via pol3 and try to correct
    Double_t y[1000], x[1000];
    Int_t n = 0;
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = 0;
    }
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    TCanvas* cfit_smiley = new TCanvas("cfit_smiley", "fit smiley", 100, 100, 800, 800);
    cfit_smiley->Clear();
    cfit_smiley->Divide(1, 4);
    cfit_smiley->cd(1);
    TH2F* histo1 = (TH2F*)histo->Clone();
    histo1->Draw("colz");
    TH2F* histo2 = (TH2F*)histo->Clone();
    histo2->RebinX(50);
    histo2->GetYaxis()->SetRangeUser(fTofdToTLow, fTofdToTHigh);
    // histo2->SetAxisRange(fTofdTotLow,fTofdTotHigh,"Y");
    cfit_smiley->cd(2);
    histo2->Draw("colz");
    std::cout << "Searching for points to fit...\n";
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        // std::cout<<"Bin "<<i<<" of "<<histo2->GetNbinsX()<<" with cut: "<<fTofdTotLow<<" < sqrt(q1*q2) <
        // "<<fTofdTotHigh<<"\n";
        cfit_smiley->cd(2);
        TLine* l = new TLine(
            histo2->GetXaxis()->GetBinCenter(i), fTofdToTLow, histo2->GetXaxis()->GetBinCenter(i), fTofdToTHigh);
        l->SetLineColor(kRed);
        l->SetLineWidth(2.);
        l->Draw();
        cfit_smiley->cd(3);
        TH1F* histo_py = (TH1F*)histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        // cfit_smiley->Update();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);

        if ((x[n] < min || x[n] > max) || (y[n] < fTofdToTLow || y[n] > fTofdToTHigh))
        {
            delete histo_py;
            continue;
        }
        if (histo_py->GetMaximum() > 5)
        {
            n++;
            delete l;
        }
        delete histo_py;
    }
    gr1 = new TGraph(n, x, y);
    gr1->SetTitle("Points found for fitting; x position in cm; sqrt(tot1*tot2)");
    gr1->Draw("A*");
    std::cout << "Start fitting\n";
    TF1* f1 = new TF1("f1", "pol3", min, max);
    f1->SetLineColor(2);
    gr1->Fit("f1", "Q", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f1->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    // fit again but with more information and better cuts
    std::cout << "Fit again with more information\n";
    n = 0;
    cfit_smiley->cd(4);
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        Double_t pos = histo2->GetXaxis()->GetBinCenter(i);
        Double_t ymean = f1->Eval(pos);
        histo2->SetAxisRange(ymean - 5., ymean + 5., "Y");
        histo2->Draw("colz");
        TH1F* histo_py = (TH1F*)histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 2)
            n++;
        delete histo_py;
    }
    gr2 = new TGraph(n, x, y);
    gr2->SetTitle("More information;x position in cm;sqrt(q1*q2)");
    gr2->Draw("A*");
    f1->DrawCopy("SAME");
    TF1* f2 = new TF1("f2", "pol3", min, max);
    f2->SetParameters(para[0], para[1], para[2], para[3]);
    f2->SetLineColor(3);
    gr2->Fit("f2", "0Q", "", min, max);
    f2->Draw("SAME");
    std::cout << "Will write:\n";
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f2->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    histo2->GetYaxis()->SetRangeUser(fTofdToTLow, fTofdToTHigh);
    auto legend = new TLegend(.9, 0.7, .99, 0.9);
    legend->AddEntry("f1", "First Fit", "l");
    legend->AddEntry("f2", "Second Fit", "l");
    legend->Draw();
    cfit_smiley->Update();
    // gPad->WaitPrimitive();
    gSystem->Sleep(3000);
    delete histo1;
    delete histo2;
    delete gr1;
    delete gr2;
    delete f1;
    delete f2;
    delete cfit_smiley;
}

void R3BTofdCal2HitParS473::fitWalk(TH2F* histo, Double_t min, Double_t max, Double_t* para){
   Double_t x[1000], y[1000];
   Int_t n = 0;
   for(Int_t j = 0; j < 5; j++){
   	para[j] = 0.;
   }
   TGraph* gr1 = new TGraph();
   TGraph* gr2 = new TGraph();
   TCanvas* cFitWalk = new TCanvas("cfit_Walk","fit walk",100,100,800,800);
   cFitWalk->Clear();
   cFitWalk->Divide(1,3);
   cFitWalk->cd(1);
   TH2F* histo1 = (TH2F*) histo->Clone();
   TH2F* histo2 = (TH2F*) histo->Clone();
   histo1->Draw("colz");
   cFitWalk->cd(2); 
}

void R3BTofdCal2HitParS473::setWalkPar()
{
	for(Int_t i = 0; i < fNofPlanes; i++){
		for(Int_t j = 0; j < fPaddlesPerPlane; j++){
			R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
			par->SetPar1Walk(5.121991e+01);
			par->SetPar2Walk(-9.9738836e-02);
			par->SetPar3Walk(5.817890e+01);
			par->SetPar4Walk(8.208630e-03);
			par->SetPar5Walk(-2.719159e-05);
		}
	}

	fCal_Par->setChanged();
}

Double_t R3BTofdCal2HitParS473::walk(Double_t Q)
{
    Double_t y = 0;
    Double_t par1, par2, par3, par4, par5;

    Int_t voltage = 443;

    if (voltage == 443)
    {
        par1 = 5.121991e+01;
        par2 = -9.9738836e-02;
        par3 = 5.817890e+01;
        par4 = 8.208630e-03;
        par5 = -2.719159e-05;
    }

    if (voltage == 444)
    {
        par1 = 1.16631e+01;
        par2 = 3.03014e-01;
        par3 = 3.25826e+02;
        par4 = -2.04208e-01;
        par5 = 3.47279e-04;
    }

    if (voltage == 500)
    {
        par1 = 1.64344e+01;
        par2 = 2.84000e-01;
        par3 = 3.47659e+02;
        par4 = -2.70050e-01;
        par5 = 3.61515e-04;
    }

    if (voltage == 600)
    {
        par1 = 1.22606e+01;
        par2 = 3.12697e-01;
        par3 = 4.40109e+02;
        par4 = -1.86328e-01;
        par5 = 1.49519e-04;
    }

    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;

    return y;
}

Double_t R3BTofdCal2HitParS473::saturation(Double_t x)
{
    Double_t kor;
    Int_t voltage = 600;

    if (voltage == 600)
    {
        if (x < 173)
        {
            kor = 0.;
        }
        else if (x > 208)
        {
            kor = -1.73665e+03 + 2.82009e+01 * 208. - 1.53846e-01 * (208. * 208.) + 2.82425e-04 * (208. * 208. * 208.);
        }
        else
        {
            kor = -1.73665e+03 + 2.82009e+01 * x - 1.53846e-01 * (x * x) + 2.82425e-04 * (x * x * x);
        }
    }

    if (voltage == 500)
    {
        if (x < 95.5)
        {
            kor = 0.;
        }
        else if (x > 124)
        {
            kor = 1.08 * x - 112.44;
        }
        else
        {
            kor = 643.257 - 16.7823 * x + 0.139822 * (x * x) - 0.000362154 * (x * x * x);
        }
    }

    if (voltage == 700)
    {
        if (x < 198)
        {
            kor = 0.;
        }
        else if (x > 298)
        {
            kor = 0.21 * x - 45.54;
        }
        else
        {
            kor = -19067 + 383.93 * x - 3.05794 * (x * x) + 0.0120429 * (x * x * x) - 2.34619e-05 * (x * x * x * x) +
                  1.81076e-08 * (x * x * x * x * x);
        }
    }

    return kor;
}

ClassImp(R3BTofdCal2HitParS473)
