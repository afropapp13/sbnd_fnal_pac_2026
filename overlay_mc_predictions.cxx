#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "constants.h"
#include "helper_functions.h"

using namespace std;
using namespace constants;

//----------------------------------------//

void overlay_mc_predictions() {

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		
	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.07,"t");
	gStyle->SetOptStat(0);

	//----------------------------------------//

	// AR23 uncertainty file

	/*TFile* ar23_unc_file = TFile::Open("mc_files/covariances.root","readonly");*/

	//----------------------------------------//

	// grab one file and create list of plots
	TFile* f = TFile::Open("output_files/analyzer_AR23.root","readonly");

	vector<TString> PlotNames = get_th1d_names(f); 

	const int nplots = PlotNames.size();
	cout << "Number of 1D Plots = " << nplots << endl;

	//----------------------------------------//

	// 1st index = histo name, 2nd index = sample
	vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear(); PlotsTrue.resize(nplots);						

	vector<TString> NameOfSamples; NameOfSamples.clear();
	vector<int> Colors; Colors.clear();		
	vector<TString> Labels; Labels.clear();
	vector<int> LineStyle; LineStyle.clear();
	
	//----------------------------------------//	

	NameOfSamples.push_back("AR23"); Colors.push_back(kOrange+7); Labels.push_back("AR23 hA");LineStyle.push_back(kSolid);
	NameOfSamples.push_back("AR25"); Colors.push_back(kGreen+2); Labels.push_back("AR25 hA"); LineStyle.push_back(kSolid);
	NameOfSamples.push_back("AR25_20j"); Colors.push_back(kBlue+1); Labels.push_back("AR25 hN");LineStyle.push_back(kSolid);
	NameOfSamples.push_back("AR25_20l"); Colors.push_back(kMagenta+1); Labels.push_back("AR25 G4");LineStyle.push_back(kSolid);
	NameOfSamples.push_back("NEUT"); Colors.push_back(kAzure+7); Labels.push_back("NEUT");LineStyle.push_back(kSolid);			

	//----------------------------------------//

	const int nsamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();
	vector<TH2D*> cov_matrix; cov_matrix.clear(); cov_matrix.resize(nplots);

	//----------------------------------------//

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot ++) {	

		//----------------------------------------//		

		TCanvas* PlotCanvas = new TCanvas(PlotNames[iplot],PlotNames[iplot],205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetBottomMargin(0.12);
		PlotCanvas->SetTopMargin(0.09);
		PlotCanvas->SetLeftMargin(0.16);
		PlotCanvas->SetRightMargin(0.05);	
		PlotCanvas->cd();			

		TLegend* leg = new TLegend(0.37,0.75,0.92,0.89);
		leg->SetBorderSize(0);
		leg->SetTextSize(text_size);
		leg->SetTextFont(text_font);
		leg->SetNColumns(2);
		leg->SetMargin(0.15);			
		leg->SetFillStyle(0);   // transparent background
		leg->SetFillColor(0);   // just in case			

		//----------------------------------------//

		PlotsTrue.at(iplot).resize(nsamples);

		for (int isample = 0; isample < nsamples; isample ++) {	

			TFile* f = TFile::Open("output_files/analyzer_"+NameOfSamples[isample]+".root");											
			PlotsTrue[iplot][isample] = (TH1D*)(f->Get(PlotNames[iplot]));	

			PlotsTrue[iplot][isample]->GetXaxis()->CenterTitle();
			PlotsTrue[iplot][isample]->GetYaxis()->CenterTitle();
			
			PlotsTrue[iplot][isample]->GetXaxis()->SetTitleFont(text_font);
			PlotsTrue[iplot][isample]->GetYaxis()->SetTitleFont(text_font);

			PlotsTrue[iplot][isample]->GetYaxis()->SetTitleOffset(1.2);

			gStyle->SetTitleFont(text_font,"");
			gStyle->SetTitleSize(text_size,"");			

			PlotsTrue[iplot][isample]->SetLineColor(Colors[isample]);
			PlotsTrue[iplot][isample]->SetLineStyle(LineStyle[isample]);
			PlotsTrue[iplot][isample]->SetMarkerColor(Colors[isample]);
			PlotsTrue[iplot][isample]->SetLineWidth(3);
			
			double max = TMath::Max(PlotsTrue[iplot][isample]->GetMaximum(), PlotsTrue[iplot][0]->GetMaximum());
			PlotsTrue[iplot][0]->SetMaximum(1.1*max);

			PlotsTrue[iplot][0]->Draw("hist same");			
			PlotsTrue[iplot][isample]->Draw("hist same");	

			double Chi2[nsamples];
			int Ndof[nsamples];
			double pval[nsamples];		
			
		//----------------------------------------//

		// only for AR23, add systematics

		if (NameOfSamples[isample] == "AR23") {

// this needs to go away when the correct covariance matrix is implemented for all the samples, including AR23			
int n = PlotsTrue[iplot][isample]->GetNbinsX();
TH2D* hCov = new TH2D("hCov","Covariance",
                      n,0,n,
                      n,0,n);

for(int i=1;i<=n;i++){

    for(int j=1;j<=n;j++){

		if (i != j) {hCov->SetBinContent(i,j,0.); }
		else {hCov->SetBinContent(i,j,0.01*PlotsTrue[iplot][0]->GetBinContent(j)); }

    }
}	

cov_matrix[iplot] = hCov;

				//this needs to be added back
				/*int n = PlotsTrue[iplot][isample]->GetNbinsX();
				TString cov_name = "tot_covariance_"+PlotNames[iplot];
				cov_matrix[iplot] = (TH2D*)ar23_unc_file->Get(cov_name);	
				cov_matrix[iplot]->Scale(1e-80);	
			
				divide_bin_area(cov_matrix[iplot],PlotsTrue[iplot][isample]);*/
				set_unc_from_cov(cov_matrix[iplot],PlotsTrue[iplot][0]);

				//----------------------//

				// error band on AR23 prediction

				TH1D* h = PlotsTrue[iplot][isample];

				// Build graph with symmetric errors
				TGraphAsymmErrors* gerr = new TGraphAsymmErrors(h);

				for (int i = 0; i < gerr->GetN(); ++i) {
					double x, y;
					gerr->GetPoint(i, x, y);

					// True symmetric error
					double eyup  = gerr->GetErrorYhigh(i);
					double eydn  = gerr->GetErrorYlow(i);
					double ey    = 0.5 * (eyup + eydn);  // enforce symmetry

					// Keep point at central value
					gerr->SetPoint(i, x, y);

					// Make symmetric around y
					gerr->SetPointEYhigh(i, ey);
					gerr->SetPointEYlow(i, ey);
				}

				// Style only the symmetric band
				gerr->SetFillColor(Colors[isample]);
				gerr->SetFillStyle(3004);
				gerr->SetLineWidth(0);

				// Draw histogram line
				h->SetLineColor(Colors[isample]);
				h->SetFillStyle(0);
				h->Draw("hist same");

				// Draw symmetric error band ONLY
				gerr->Draw("E2 same");

		} // end of AR23 

			calc_chi2(PlotsTrue[iplot][0],PlotsTrue[iplot][isample],cov_matrix[iplot],Chi2[isample],Ndof[isample],pval[isample]);
			TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[isample],1) + "/" + TString(std::to_string(Ndof[isample])) +")";

			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[iplot][isample],Labels[isample] + Chi2NdofAlt,"");
			lGenie->SetTextColor(Colors[isample]); 

		} // end of loop over samples

		//----------------------------------------//

		TLatex tex;
		tex.SetTextFont(text_font);
		tex.SetTextSize(text_size);	
		tex.SetTextColor(kGray+1);					
		tex.DrawLatexNDC(0.18, 0.86, "#bf{SBND} Internal");	

		leg->Draw("same");
		PlotCanvas->SaveAs("plots/mc_prediction_"+PlotNames[iplot]+".pdf");
		delete PlotCanvas;

		//----------------------------------------//		

	} // End of the loop over the plots

} // End of the program 
