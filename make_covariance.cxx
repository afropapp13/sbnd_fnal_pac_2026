#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include "helper_functions.h"

using namespace std;

//--------------------//

void make_covariance() {

	//--------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptStat(0);

	//--------------------//

	vector< tuple<int,int>> sigma_map = {
		{0, 1},
		{1, -1},
		{2, 2},
		{3, -2},
		{4, 3},								
		{5, -3},
		{6, 0}
	};	

	//--------------------//

	// cv mc sample

	TFile* cv_mc_file = TFile::Open("output_files/analyzer_AR23.root","readonly");

	vector<TString> plot_names = get_th1d_names(cv_mc_file); 
	const int nplots = plot_names.size();
	cout << "Number of 1D Plots = " << nplots << endl;	

	// cv plots

	vector <TH1D*> cv_plots; cv_plots.resize(nplots);

	for (int i = 0; i < nplots; i++) {

		cv_plots[i] = (TH1D*)(cv_mc_file->Get(plot_names[i]));
		cv_plots[i]->SetDirectory(0); // to decouple it from the open file directory

	}

	cv_mc_file->Close();

	//--------------------//
	
	// declaration of matrices & canvases
	//1st index = plot, 2nd index = knob

	vector<TH2D*> tot_covariance; tot_covariance.resize(nplots);	
	vector< vector<TH2D*> > covariance; covariance.resize(nplots);
	vector< vector<TCanvas*> > canvas; canvas.resize(nplots);	
	vector< vector<TLegend*> > leg; leg.resize(nplots);			

	//--------------------//

	// color palette

	std::vector<Color_t> colors = {kBlue+2, kOrange+7, kGreen+2, kRed+2, kAzure+7, kMagenta+1, kGray+1};	

	//--------------------//

	// alternative mc

	std::vector<int> universe;
	std::vector<TString> knob;

	// ---------------- ZExp PCA (7 universes) ----------------
	/*knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b1"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b2"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b3"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b4"); universe.push_back(7);*/

	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b2"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b3"); universe.push_back(7);
	knob.push_back("ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b4"); universe.push_back(7);

	// ---------------- CCQE Template (SF) ----------------
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin1"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin2"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin3"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin4"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin5"); universe.push_back(6);

	// ---------------- CCQE Template (HF) ----------------
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin1"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin2"); universe.push_back(6);
	/*knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin3"); universe.push_back(6);*/
	/*knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin4"); universe.push_back(6);*/
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin5"); universe.push_back(6);

	// ---------------- CCQE Template (CRPA) ----------------
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin1"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin2"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin3"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin4"); universe.push_back(6);
	knob.push_back("CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin5"); universe.push_back(6);

	// ---------------- QE Interference ----------------
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_0"); universe.push_back(6);
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_1"); universe.push_back(6);
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_2"); universe.push_back(6);
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_3"); universe.push_back(6);
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_4"); universe.push_back(6);
	knob.push_back("QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_5"); universe.push_back(6);

	// ---------------- GENIE EDepFSI knobs ----------------
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_VecFFCCQEshape"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_CoulombCCQE"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormCCMEC"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormNCMEC"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_DecayAngMEC"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFP_pi"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrCEx_pi"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrInel_pi"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrAbs_pi"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrPiProd_pi"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCL_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4LoE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLLoE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M1E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM1E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M2E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM2E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4HiE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLHiE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPLoE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM1E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM2E_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPHiE_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProFix_N"); universe.push_back(6);
	knob.push_back("GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProBias_N"); universe.push_back(6);

	// ---------------- MEC Valencia ----------------
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin0"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin1"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin2"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin3"); universe.push_back(6);

	// ---------------- MEC Martini ----------------
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin0"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin1"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin2"); universe.push_back(6);
	knob.push_back("MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin3"); universe.push_back(6);

	//--------------------//	

	// alternative mc files and plots

	int nknobs = knob.size();

	// 1st index = knob, 2nd index = universe
	vector<TFile*> knob_file; knob_file.resize(nknobs);

	//1st index = plot, 2nd index = knob, 3rd index = universe
	vector< vector <vector <TH1D*> > > knob_plots; knob_plots.resize(nplots);

	// loop over the plots
	for (int iplot = 0; iplot < nplots; iplot++) {

		knob_plots[iplot].resize(nknobs);

		// loop over the knobs
		for (int alt = 0; alt < nknobs; alt++ ) {

			knob_plots[iplot][alt].resize( universe.at(alt) );

			// loop over the universes
			for (int iuni = 0; iuni < universe.at(alt); iuni++ ) {		

				TString TStringAltMC = "output_files/"+knob[alt]+"_" + int_to_string(iuni)+"_analyzer_AR23.root";
				knob_file[alt] = TFile::Open(TStringAltMC,"readonly");			

				knob_plots[iplot][alt][iuni] = (TH1D*)(knob_file[alt]->Get(plot_names[iplot]));
				knob_plots[iplot][alt][iuni]->SetDirectory(0); // to decouple it from the open file directory
				knob_file[alt]->Close();

			} // end of the loop over the universes

		} // end of the loop over the knobs

	} // end of the loop over the plots

	//--------------------//

	// file to store all the covariances

	TFile* cov_file = new TFile("output_files/covariances.root","recreate");	

	//--------------------//	

	// looping over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {

		canvas.at(iplot).resize(nknobs);
		leg.at(iplot).resize(nknobs);	
		covariance.at(iplot).resize(nknobs);	
		
		tot_covariance.at(iplot) = new TH2D("tot_covariance_"+plot_names[iplot], "tot_covariance_"+plot_names[iplot],
									   cv_plots[iplot]->GetNbinsX(), 0.5, cv_plots[iplot]->GetNbinsX()+0.5,
									   cv_plots[iplot]->GetNbinsX(), 0.5, cv_plots[iplot]->GetNbinsX()+0.5);		

		tot_covariance.at(iplot)->Reset();  // sets all bin contents (and errors) to 0						   

		// looping over the knobs

		for (int iknob = 0; iknob < nknobs; iknob++) {

			covariance.at(iplot).at(iknob) = new TH2D("covariance_"+knob[iknob]+"_"+plot_names[iplot], "covariance_"+knob[iknob]+"_"+plot_names[iplot],
									   cv_plots[iplot]->GetNbinsX(), 0.5, cv_plots[iplot]->GetNbinsX()+0.5,
									   cv_plots[iplot]->GetNbinsX(), 0.5, cv_plots[iplot]->GetNbinsX()+0.5);

			covariance.at(iplot).at(iknob)->Reset();  // sets all bin contents (and errors) to 0						   

			// create canvas
			canvas.at(iplot).at(iknob) = new TCanvas("c_"+knob[iknob]+"_"+plot_names[iplot],"c_"+knob[iknob]+"_"+plot_names[iplot],205,34,1024,768);
			canvas.at(iplot).at(iknob)->cd();
			canvas.at(iplot).at(iknob)->SetBottomMargin(0.15);
			canvas.at(iplot).at(iknob)->SetLeftMargin(0.15);
			canvas.at(iplot).at(iknob)->SetRightMargin(0.1);	
			canvas.at(iplot).at(iknob)->SetTopMargin(0.11);

			leg.at(iplot).at(iknob) = new TLegend(0.15,0.9,0.9,0.98);
			leg.at(iplot).at(iknob)->SetBorderSize(0);
			leg.at(iplot).at(iknob)->SetTextFont(text_font);
			leg.at(iplot).at(iknob)->SetTextSize(text_size);
			leg.at(iplot).at(iknob)->SetNColumns(7);			

			// plot cv

			cv_plots[iplot]->GetXaxis()->CenterTitle();
			cv_plots[iplot]->GetXaxis()->SetTitleFont(text_font);
			cv_plots[iplot]->GetXaxis()->SetTitleSize(text_size);
			cv_plots[iplot]->GetXaxis()->SetLabelFont(text_font);
			cv_plots[iplot]->GetXaxis()->SetLabelSize(text_size);
			cv_plots[iplot]->GetXaxis()->SetNdivisions(6);

			cv_plots[iplot]->GetYaxis()->CenterTitle();
			cv_plots[iplot]->GetYaxis()->SetTitleFont(text_font);
			cv_plots[iplot]->GetYaxis()->SetTitleSize(text_size);
			cv_plots[iplot]->GetYaxis()->SetLabelFont(text_font);
			cv_plots[iplot]->GetYaxis()->SetLabelSize(text_size);
			cv_plots[iplot]->GetYaxis()->SetNdivisions(6);
			cv_plots[iplot]->GetYaxis()->SetTitle("cross section [no bin width division]");			

			cv_plots[iplot]->SetLineWidth(3);			
			cv_plots[iplot]->SetLineColor(kBlack);
			cv_plots[iplot]->Draw("hist");
			leg.at(iplot).at(iknob)->AddEntry(cv_plots[iplot], "CV", "");		
			// loop over the sigma universes

			//for (int iuni = 0; iuni < 1; iuni++ ) {	
			for (int iuni = 0; iuni < universe.at(iknob); iuni++ ) {		

				knob_plots[iplot][iknob][iuni]->SetLineWidth(3);				
				knob_plots[iplot][iknob][iuni]->SetLineColor( colors.at(iuni) );
				knob_plots[iplot][iknob][iuni]->Draw("hist same");

				// identify the variation in sigmas and add to legend
				int sigmas = std::get<1>( sigma_map.at(iuni) );
				TString leg_label = int_to_string( sigmas ) + "#sigma";

				TLegendEntry* legColor = leg.at(iplot).at(iknob)->AddEntry(knob_plots[iplot][iknob][iuni], leg_label, "");
				legColor->SetTextColor( colors.at(iuni) ); 				

				// construct the covariance matrix for the +1 variation
				if (sigmas == 1) {

					// loop over the x bins
					for (int i = 1; i <= cv_plots[iplot]->GetNbinsX(); i++) {

						for (int j = 1; j <= cv_plots[iplot]->GetNbinsX(); j++) {

							double diff_i = knob_plots[iplot][iknob][iuni]->GetBinContent(i) - cv_plots[iplot]->GetBinContent(i);
							double diff_j = knob_plots[iplot][iknob][iuni]->GetBinContent(j) - cv_plots[iplot]->GetBinContent(j);

							double cov_ij = diff_i * diff_j;

							covariance.at(iplot).at(iknob)->SetBinContent(i, j, cov_ij);
							
						}

					}

					// storing the covariance matrix
					cov_file->cd();
					covariance.at(iplot).at(iknob)->Write("covariance_"+knob[iknob]+"_"+plot_names[iplot]);
					tot_covariance.at(iplot)->Add( covariance.at(iplot).at(iknob) );			

				}

			} // end of the loop over the sigma universes

			double max_y = cv_plots[iplot]->GetMaximum() * 1.04;
			cv_plots[iplot]->SetMaximum( max_y );
			cv_plots[iplot]->SetMinimum( 0 );			

			cv_plots[iplot]->Draw("hist same");			

			leg.at(iplot).at(iknob)->Draw("same");

			TLatex tex;
			tex.SetTextFont(text_font);
			tex.SetTextSize(text_size);			
			tex.DrawLatexNDC(0.18, 0.83, knob[iknob]);

			// add the uncertainty on the CV

			TH1D* h = cv_plots[iplot];

			// Build graph with symmetric errors
			TGraphAsymmErrors* gerr = new TGraphAsymmErrors(h);

			for (int i = 0; i < gerr->GetN(); ++i) {

				double x, y;
				gerr->GetPoint(i, x, y);

				// True symmetric error
				// double eyup  = gerr->GetErrorYhigh(i);
				// double eydn  = gerr->GetErrorYlow(i);
				// double ey    = 0.5 * (eyup + eydn);  // enforce symmetry

				// Keep point at central value
				gerr->SetPoint(i, x, y);

				double unc = TMath::Sqrt( covariance.at(iplot).at(iknob)->GetBinContent(i+1, i+1) );

				// Make symmetric around y
				gerr->SetPointEYhigh(i, unc);
				gerr->SetPointEYlow(i, unc);
			}

			// Style only the symmetric band
			gerr->SetFillColor(kRed);
			gerr->SetFillStyle(3004);
			gerr->SetLineWidth(0);

			// Draw histogram line
			h->SetLineColor(kRed);
			h->SetFillStyle(0);
			h->Draw("hist same");			

			// Draw symmetric error band ONLY
			gerr->Draw("E2 same");

			canvas.at(iplot).at(iknob)->SaveAs("pdf/knob_"+knob[iknob]+"_"+plot_names[iplot]+".pdf");

		} // end of the loop over the knobs

		cov_file->cd();
		tot_covariance.at(iplot)->Write("tot_covariance_"+plot_names[iplot]);
		cv_plots[iplot]->Write("cv_plot_"+plot_names[iplot]);

	} // end of the loop over the plots

	//--------------------//	

} // end of the program
