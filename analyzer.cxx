#define analyzer_cxx
#include "analyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <vector>
#include <cmath>

#include "constants.h"
#include "STV_Tools.h"
#include "Tools.h"
#include "helper_functions.h"

using namespace std;
using namespace constants;

void analyzer::Loop() {

    //--------------------------------------//    

    TH1::SetDefaultSumw2();
    if (!fChain) return;
    const Long64_t nentries = fChain->GetEntriesFast();

    Tools tools;
    force_plot_style();

    //--------------------------------------//
    
    TString outname = "output_files/analyzer_" + fOutputFile + ".root";

    if (fOutputFile == "AR23" && fweights != "" && findex != -1)
        outname = "output_files/" + fweights + "_" +
                  int_to_string(findex) +
                  "_analyzer_" + fOutputFile + ".root";

    TFile* outfile = TFile::Open(outname, "RECREATE");

    //--------------------------------------
  
    TFile* syst_file = nullptr;
    TTree* syst_tree = nullptr;

    double tweak_responses[7] = {1,1,1,1,1,1,1};
    int ntweaks = 0;
    double paramCVWeight;    

    if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

        syst_file = TFile::Open(
            "/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/sbnd_fnal_pac_2026/AR23_20i_00_000//syst_AR23_20i_00_000.root",
            "READonly");

        syst_tree = (TTree*)syst_file->Get("events");

        syst_tree->SetBranchStatus("*", 0);
        syst_tree->SetBranchStatus(("tweak_responses_" + fweights).Data(), 1);
        syst_tree->SetBranchStatus(("ntweaks_" + fweights).Data(), 1);
        syst_tree->SetBranchStatus(("paramCVWeight_" + fweights).Data(), 1);         

        syst_tree->SetBranchAddress(("ntweaks_" + fweights).Data(), &ntweaks);
        syst_tree->SetBranchAddress(("tweak_responses_" + fweights).Data(), &tweak_responses);
        syst_tree->SetBranchAddress(("paramCVWeight_" + fweights).Data(), &paramCVWeight);

    }

    //--------------------------------------//

    outfile->cd();
    
    // 1d plots

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
        ";#deltap_{T};#frac{d#sigma}{d#deltap_{T}} [10^{-38} cm^{2}/(Ar GeV/c)]",
        NBinsDeltaPT,
        ArrayNBinsDeltaPT
    );

    TH1D* DeltaAlphaTPlot = new TH1D(
        "DeltaAlphaTPlot",
        ";#delta#alpha_{T};#frac{d#sigma}{d#delta#alpha_{T}} [10^{-38} cm^{2}/(Ar deg)]",
        NBinsDeltaAlphaT,
        ArrayNBinsDeltaAlphaT
    ); 

    TH1D* ProtonMultiPlot = new TH1D(
        "ProtonMultiPlot",
        ";n_{p};#frac{d#sigma}{dn_{p}} [10^{-38} cm^{2}/(Ar)]",
        nproton_bins,ArrayNBinsProtonMulti);  

    TH1D* PionMultiPlot = new TH1D(
        "PionMultiPlot",
        ";n_{#pi};#frac{d#sigma}{dn_{#pi}} [10^{-38} cm^{2}/(Ar)]",
        npion_bins,ArrayNBinsPionMulti);        

    //--------------------------------------//

    // 3d uncorrelated plots

    TH1D* dpt_in_dat_proton_multiPlot[nproton_bins][TwoDNBinsDeltaAlphaT];
    TH1D* pmu_in_costheta_mu_proton_multiPlot[nproton_bins][TwoDNBinsMuonCosTheta];
    TH1D* pmu_in_costheta_mu_pion_multiPlot[npion_bins][TwoDNBinsMuonCosTheta];        

    // loop over proton multiplicity
    for (int i_prot = 0; i_prot < (int)nproton_bins; ++i_prot) {

        // dpt in slices of dat and proton multiplicity
        for (int i_dat = 0; i_dat < TwoDNBinsDeltaAlphaT; ++i_dat) {

            TString name = "DeltaPT_in_ProtonMulti_" + int_to_string(i_prot) + "_DeltaAlphaT_" + int_to_string(i_dat);
            TString title = int_to_string(i_prot) + " protons, " + to_string_with_precision(TwoDArrayNBinsDeltaAlphaT[i_dat],0) + " < #delta#alpha_{T} < " + to_string_with_precision(TwoDArrayNBinsDeltaAlphaT[i_dat+1],0) + " deg;#deltap_{T} [GeV/c];#frac{d^{2}#sigma}{d#deltap_{T}d#delta#alpha_{T}dn_{p}} [10^{-38} cm^{2}/(Ar GeV/c deg)]";
            dpt_in_dat_proton_multiPlot[i_prot][i_dat] = new TH1D(name, title, bins_dpt_in_dat_proton_mult.at(i_prot).at(i_dat).size() - 1, bins_dpt_in_dat_proton_mult.at(i_prot).at(i_dat).data());

        } // end of dpt in slices of dat and proton multiplicity

        // pmu in slices of costheta_mu and proton multiplicity
        for (int i_costheta_mu = 0; i_costheta_mu < TwoDNBinsMuonCosTheta; ++i_costheta_mu) {

            TString name = "Pmu_in_ProtonMulti_" + int_to_string(i_prot) + "_MuonCosTheta_" + int_to_string(i_costheta_mu);
            TString title = int_to_string(i_prot) + " protons, " + to_string_with_precision(TwoDArrayNBinsMuonCosTheta[i_costheta_mu],2) + " < cos(#theta_{#mu}) < " + to_string_with_precision(TwoDArrayNBinsMuonCosTheta[i_costheta_mu+1],2) + ";p_{#mu} [GeV/c];#frac{d^{2}#sigma}{dp_{#mu}dcos(#theta_{#mu})dn_{p}} [10^{-38} cm^{2}/(Ar GeV/c)]";
            pmu_in_costheta_mu_proton_multiPlot[i_prot][i_costheta_mu] = new TH1D(name, title, bins_pmu_in_mucostheta_proton_mult.at(i_prot).at(i_costheta_mu).size() - 1, bins_pmu_in_mucostheta_proton_mult.at(i_prot).at(i_costheta_mu).data());

        } // end of dpt in slices of dat and proton multiplicity        

    } // end of the loop over proton multiplicity

    // loop over pion multiplicity
    for (int i_pi = 0; i_pi < (int)npion_bins; ++i_pi) {    

        // pmu in slices of costheta_mu and pion multiplicity
        for (int i_costheta_mu = 0; i_costheta_mu < TwoDNBinsMuonCosTheta; ++i_costheta_mu) {

            TString name = "Pmu_in_PionMulti_" + int_to_string(i_pi) + "_MuonCosTheta_" + int_to_string(i_costheta_mu);
            TString title = int_to_string(i_pi) + " pions, " + to_string_with_precision(TwoDArrayNBinsMuonCosTheta[i_costheta_mu],2) + " < cos(#theta_{#mu}) < " + to_string_with_precision(TwoDArrayNBinsMuonCosTheta[i_costheta_mu+1],2) + ";p_{#mu} [GeV/c];#frac{d^{2}#sigma}{dp_{#mu}dcos(#theta_{#mu})dn_{#pi}} [10^{-38} cm^{2}/(Ar GeV/c)]";
            pmu_in_costheta_mu_pion_multiPlot[i_pi][i_costheta_mu] = new TH1D(name, title, bins_pmu_in_mucostheta_pion_mult.at(i_pi).at(i_costheta_mu).size() - 1, bins_pmu_in_mucostheta_pion_mult.at(i_pi).at(i_costheta_mu).data());

        } // end of pmu in slices of costheta_mu and pion multiplicity        

    } // end of the loop over pion multiplicity

    //--------------------------------------//
    
    // 3d correlated plots

    TH1D* serial_dpt_in_dat_proton_multiPlot  = \
        new TH1D("serial_dpt_in_dat_proton_multiPlot",\
        xlabel_serial_dpt_in_dat_proton_multi,\
        tools.Return3DNBins(bins_dpt_in_dat_proton_mult),\
        &tools.Return3DBinIndices(bins_dpt_in_dat_proton_mult)[0]);
        
    TH1D* serial_pmu_in_costheta_mu_proton_multiPlot  = \
        new TH1D("serial_pmu_in_costheta_mu_proton_multiPlot",\
        xlabel_serial_pmu_in_costhetamu_proton_multi,\
        tools.Return3DNBins(bins_pmu_in_mucostheta_proton_mult),\
        &tools.Return3DBinIndices(bins_pmu_in_mucostheta_proton_mult)[0]);     
        
    TH1D* serial_pmu_in_costheta_mu_pion_multiPlot  = \
        new TH1D("serial_pmu_in_costheta_mu_pion_multiPlot",\
        xlabel_serial_pmu_in_costhetamu_pion_multi,\
        tools.Return3DNBins(bins_pmu_in_mucostheta_pion_mult),\
        &tools.Return3DBinIndices(bins_pmu_in_mucostheta_pion_mult)[0]);             

    //--------------------------------------//

    const int reportEvery = std::max((Long64_t)1, nentries / 100);   
    cout << endl << outname << " processing" << endl; 
    
    for (Long64_t i = 0; i < nentries; ++i) {

        //--------------------------------------//

        // progress bar
        if (i % reportEvery == 0 || i == nentries - 1) {
            double percent = 100.0 * i / nentries;
            std::cout << "\rprocessing: "
                    << std::fixed << std::setprecision(1)
                    << percent << "% "
                    << std::flush;
        }        

        //--------------------------------------//        

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

        //----------------------------------//        

        double syst_weight = 1.0;

        if (findex >= 0) {
            
            syst_tree->GetEntry(i);    
            syst_weight = tweak_responses[findex];
            // tweak_responses[6] has the new CV
            if (fweights.Contains("MvA")) { syst_weight = syst_weight/tweak_responses[6]; }

        }

        double weight = fScaleFactor * Units * A * Weight * syst_weight;
        if (std::isnan(weight)) { continue; }

        //----------------------------------//            

        // NEUT 6.1.2 was produced in 14 batches
        if (fOutputFile == "NEUT") { weight /= 14.; }
        if (fOutputFile == "AR23" || fOutputFile == "AR25" || fOutputFile == "AR25_20j" || fOutputFile == "AR25_20l") { weight /= 10.; }         

        //----------------------------------//        

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
        double dat = stv.ReturnDeltaAlphaT();

        //--------------------------------------//

        // overflow / underflow

        if (dpt > ArrayNBinsDeltaPT[NBinsDeltaPT]) {
            dpt = 0.9999 * ArrayNBinsDeltaPT[NBinsDeltaPT];
        }

        if (npr > (int)nproton_bins-1) { npr = nproton_bins - 1; }
        if (npi > (int)npion_bins-1) { npi = npion_bins - 1; }    

        //--------------------------------------//

        // 1d plots
        MuonCosThetaPlot->Fill(mu4.Vect().CosTheta(), weight); 
        MuonMomentumPlot->Fill(mu4.Rho(), weight);
        DeltaPTPlot->Fill(dpt, weight);
        DeltaAlphaTPlot->Fill(dat, weight);
        ProtonMultiPlot->Fill(npr, weight);
        PionMultiPlot->Fill(npi, weight);

        //--------------------------------------//  
        
        // 3d uncorrelated plots

        int dat_slice = FindBin(TwoDArrayNBinsDeltaAlphaT, dat);
        int costheta_mu_slice = FindBin(TwoDArrayNBinsMuonCosTheta, mu4.Vect().CosTheta());   

        dpt_in_dat_proton_multiPlot[npr][dat_slice]->Fill(dpt, weight);
        pmu_in_costheta_mu_proton_multiPlot[npr][costheta_mu_slice]->Fill(mu4.Rho(), weight);
        pmu_in_costheta_mu_pion_multiPlot[npi][costheta_mu_slice]->Fill(mu4.Rho(), weight);

        //--------------------------------------//     
        
        // 3d correlated plots

        int dpt_index_dat_np = ReturnIndexIn3DList(bins_dpt_in_dat_proton_mult, npr, dat_slice, dpt);
        int pmu_index_cosmu_np = ReturnIndexIn3DList(bins_pmu_in_mucostheta_proton_mult, npr, costheta_mu_slice, mu4.Rho());
        int pmu_index_cosmu_npi = ReturnIndexIn3DList(bins_pmu_in_mucostheta_pion_mult, npi, costheta_mu_slice, mu4.Rho());

        serial_dpt_in_dat_proton_multiPlot->Fill(dpt_index_dat_np,weight);
        serial_pmu_in_costheta_mu_proton_multiPlot->Fill(pmu_index_cosmu_np,weight);
        serial_pmu_in_costheta_mu_pion_multiPlot->Fill(pmu_index_cosmu_npi,weight);
        
    }

    //--------------------------------------//

    // no need to divide by scaling factor 
    //for N-dim slice for 1D plots
    divide_bin_width(MuonCosThetaPlot); 
    divide_bin_width(DeltaPTPlot);
    divide_bin_width(DeltaAlphaTPlot);
    divide_bin_width(MuonMomentumPlot);

    //--------------------------------------//
    
    // 3d plots need the scaling factor

    // loop over proton multiplicity
    for (int i_prot = 0; i_prot < (int)nproton_bins; ++i_prot) {

        // dpt in slices of dat and proton multiplicity
        for (int i_dat = 0; i_dat < TwoDNBinsDeltaAlphaT; ++i_dat) {

            dpt_in_dat_proton_multiPlot[i_prot][i_dat]->Scale( 1. / (TwoDArrayNBinsDeltaAlphaT[i_dat+1] - TwoDArrayNBinsDeltaAlphaT[i_dat]) );
            divide_bin_width(dpt_in_dat_proton_multiPlot[i_prot][i_dat]);

        } // end of dpt in slices of dat and proton multiplicity

        // pmu in slices of costheta_mu and proton multiplicity
        for (int i_costheta_mu = 0; i_costheta_mu < TwoDNBinsMuonCosTheta; ++i_costheta_mu) {

            pmu_in_costheta_mu_proton_multiPlot[i_prot][i_costheta_mu]->Scale( 1. / (TwoDArrayNBinsMuonCosTheta[i_costheta_mu+1] - TwoDArrayNBinsMuonCosTheta[i_costheta_mu]) );
            divide_bin_width(pmu_in_costheta_mu_proton_multiPlot[i_prot][i_costheta_mu]);

        } // end of dpt in slices of dat and proton multiplicity        

    } // end of the loop over proton multiplicity

    // loop over pion multiplicity
    for (int i_pi = 0; i_pi < (int)npion_bins; ++i_pi) {    

        // pmu in slices of costheta_mu and pion multiplicity
        for (int i_costheta_mu = 0; i_costheta_mu < TwoDNBinsMuonCosTheta; ++i_costheta_mu) {

            pmu_in_costheta_mu_pion_multiPlot[i_pi][i_costheta_mu]->Scale( 1. / (TwoDArrayNBinsMuonCosTheta[i_costheta_mu+1] - TwoDArrayNBinsMuonCosTheta[i_costheta_mu]) );
            divide_bin_width(pmu_in_costheta_mu_pion_multiPlot[i_pi][i_costheta_mu]);

        } // end of pmu in slices of costheta_mu and pion multiplicity        

    } // end of the loop over pion multiplicity    

    //--------------------------------------//    

    outfile->Write();
    outfile->Close();
    if (syst_file) syst_file->Close();

    cout << outname << " processed" << endl;

    //--------------------------------------// 

}