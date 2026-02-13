#define analyzer_cxx
#include "analyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>
#include <vector>
#include <cmath>

#include "constants.h"
#include "STV_Tools.h"
#include "Tools.h"
#include "helper_functions.h"

using namespace std;
using namespace constants;

    //--------------------------------------//

void analyzer::Loop() {

    //--------------------------------------//    

    TH1::SetDefaultSumw2();
    if (!fChain) return;
    const Long64_t nentries = fChain->GetEntriesFast();

    Tools tools;

    //--------------------------------------//
    
    TString outname = "output_files/analyzer_" + fOutputFile + ".root";

    if (fOutputFile == "AR23" && fweights != "" && findex != -1)
        outname = "../mc_files/" + fweights + "_" +
                  int_to_string(findex) +
                  "_analyzer_" + fOutputFile + ".root";

    TFile* outfile = TFile::Open(outname, "RECREATE");

    //--------------------------------------
  
    TFile* syst_file = nullptr;
    TTree* syst_tree = nullptr;

    double tweak_responses[7] = {1,1,1,1,1,1,1};
    int ntweaks = 0;

    if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

        syst_file = TFile::Open(
            "../mc_files/syst_14_1000180400_CC_v3_6_2_AR23_20i_00_000.root",
            "READonly");

        syst_tree = (TTree*)syst_file->Get("events");

        syst_tree->SetBranchStatus("*", 0);
        syst_tree->SetBranchStatus(("tweak_responses_" + fweights).Data(), 1);
        syst_tree->SetBranchStatus(("ntweaks_" + fweights).Data(), 1);

        syst_tree->SetBranchAddress(("ntweaks_" + fweights).Data(), &ntweaks);
        syst_tree->SetBranchAddress(("tweak_responses_" + fweights).Data(), &tweak_responses);

    }

    //--------------------------------------

    outfile->cd();
    
    TH1D* MuonCosThetaPlot = new TH1D("MuonCosThetaPlot",
        ";cos(#theta_{#mu});#frac{d#sigma}{dcos(#theta_{#mu})} [10^{-38} cm^{2}/Ar]",
        NBinsMuonCosTheta,
        ArrayNBinsMuonCosTheta);

    TH1D* MuonMomentumPlot = new TH1D(
        "MuonMomentumPlot",
        ";p_{#mu} [GeV];#frac{d#sigma}{dp_{#mu}} [10^{-38} cm^{2} / (Ar GeV/c)]",
        NBinsMuonMomentum,
        ArrayNBinsMuonMomentum
    );        

    TH1D* DeltaPTPlot = new TH1D(
        "DeltaPTPlot",
        ";#deltap_{T};#frac{d#sigma}{d#deltap_{T}} [10^{-38} cm^{2}]/(Ar GeV/c)",
        NBinsDeltaPT,
        ArrayNBinsDeltaPT
    );

    TH1D* DeltaAlphaTPlot = new TH1D(
        "DeltaAlphaTPlot",
        ";#delta#alpha_{T};#frac{d#sigma}{d#delta#alpha_{T}} [10^{-38} cm^{2}/(Ar deg)]",
        NBinsDeltaAlphaT,
        ArrayNBinsDeltaAlphaT
    );

    //--------------------------------------
    
    for (Long64_t i = 0; i < nentries; ++i) {

        fChain->GetEntry(i);

        if (PDGLep != 13) continue;

        int mu = -1, maxpr = -1;
        double max_p_pr = -1.;
        int nmu = 0, npr = 0, npi = 0;

        for (int j = 0; j < nfsp; ++j) {

            const double p = std::sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j]);

            if (pdg[j] == 13) { mu = j; ++nmu; }
            else if (pdg[j] == 2212) { 

                ++npr; 
                // store the index of the most energetic proton
                if (p > max_p_pr) { max_p_pr = p; maxpr = j; }

            }
            else if (abs(pdg[j]) == 211 || pdg[j] == 111) { ++npi; }

        }

        //----------------------------------        

        double syst_weight = 1.0;

        if (findex >= 0) {
            
            syst_tree->GetEntry(i);    
            syst_weight = tweak_responses[findex];

        }

        const double weight =
            fScaleFactor * Units * A * Weight * syst_weight;

        TLorentzVector mu4(px[mu], py[mu], pz[mu], E[mu]);
        TLorentzVector pr4(0, 0, 0, 0);
        if (npr > 0) { pr4.SetPxPyPzE(px[maxpr], py[maxpr], pz[maxpr], E[maxpr]); }

        STV_Tools stv(
            mu4.Vect(),
            pr4.Vect(),
            mu4.E(),
            sqrt(pr4.Rho()*pr4.Rho() +
                 ProtonMass_GeV*ProtonMass_GeV));

        double dpt = stv.ReturnPt();
        if (dpt > ArrayNBinsDeltaPT[NBinsDeltaPT])
            dpt = 0.9999 * ArrayNBinsDeltaPT[NBinsDeltaPT];

        // 1D plots
        MuonCosThetaPlot->Fill(mu4.Vect().CosTheta(), weight); 
        MuonMomentumPlot->Fill(mu4.Rho(), weight);
        DeltaPTPlot->Fill(dpt, weight);
        DeltaAlphaTPlot->Fill(stv.ReturnDeltaAlphaT(), weight);

        
    }

    //--------------------------------------//

    // no need to divide by scaling factor for N-dim slice
    divide_bin_width(MuonCosThetaPlot); 
    divide_bin_width(DeltaPTPlot);
    divide_bin_width(DeltaAlphaTPlot);
    divide_bin_width(MuonMomentumPlot);

    //--------------------------------------//    

    outfile->Write();
    outfile->Close();
    if (syst_file) syst_file->Close();

    cout << outname << " processed" << endl;

    //--------------------------------------// 

}