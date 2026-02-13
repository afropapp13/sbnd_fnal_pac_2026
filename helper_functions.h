#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include "Math/DistFunc.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm> // Required for std::max_element

using namespace std;

//---------------------------//

std::vector<std::string> get_th1d_names(TFile* file)
{
    std::vector<std::string> names;

    if (!file || file->IsZombie()) {
        std::cerr << "Invalid file!" << std::endl;
        return names;
    }

    TIter next(file->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

        TObject* obj = key->ReadObj();

        // Check if object is a TH1D
        if (obj->InheritsFrom(TH1D::Class())) {
            names.push_back(obj->GetName());
        }
    }

    return names;
}

//---------------------------//

TString int_to_string(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;
  
}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 1) {

    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());

}

//----------------------------------------//		

void divide_bin_width(TH1D* h, double SF=1.) {

	int NBins = h->GetXaxis()->GetNbins();
  
	for (int i = 0; i < NBins; i++) {
  
	  double CurrentEntry = h->GetBinContent(i+1);
	  double NewEntry = SF * CurrentEntry / h->GetBinWidth(i+1);
  
	  double CurrentError = h->GetBinError(i+1);
	  double NewError = SF * CurrentError / h->GetBinWidth(i+1);
  
	  h->SetBinContent(i+1,NewEntry); 
	  h->SetBinError(i+1,NewError); 
	  //h->SetBinError(i+1,0.000001); 
  
	}
  
  }

//----------------------------------------//	

double Chi2Prob(double chi2, int ndof)
{
    return ROOT::Math::chisquared_cdf_c(chi2, ndof);
}

//----------------------------------------//	

void calc_chi2(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {			

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	TMatrixD inverse_cov_m = cov_m.Invert();

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;

		}

	}

	ndof = h_data_clone->GetNbinsX();
	pval = Chi2Prob(chi, ndof);

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

//----------------------------------------//	