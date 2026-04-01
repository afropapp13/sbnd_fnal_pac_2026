#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include "Math/DistFunc.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TDecompChol.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm> // Required for std::max_element
#include <cstddef>
#include "constants.h"

using namespace std;
using namespace constants;

#include <iostream>
#include <TH1D.h>

//---------------------------//

void PrintBinErrors(TH1D* h) {
    if (!h) {
        std::cout << "Histogram pointer is null!" << std::endl;
        return;
    }

    int nbins = h->GetNbinsX();

    std::cout << "Histogram: " << h->GetName() << std::endl;
    std::cout << "Bin\tCenter\t\tContent\t\tError" << std::endl;

    for (int i = 1; i <= nbins; ++i) {  // skip under/overflow for clarity
        double center  = h->GetBinCenter(i);
        double content = h->GetBinContent(i);
        double error   = h->GetBinError(i);

        std::cout << i << "\t"
                  << center << "\t"
                  << content << "\t"
                  << error << std::endl;
    }
}

//---------------------------//

void set_unc_from_cov(TH2D* cov, TH1D* h) {

	if (!cov || !h) return;

	int nbins = h->GetNbinsX();

	for (int i = 1; i <= nbins; ++i) {
		double variance = cov->GetBinContent(i, i);
		double error = (variance >= 0) ? sqrt(variance) : 0.0;
		h->SetBinError(i, error);
	}

}

//---------------------------//

void force_plot_style() {

    gStyle->SetOptStat(0);

    gStyle->SetTextFont(text_font);
    gStyle->SetLabelFont(text_font, "XYZ");
    gStyle->SetLabelFont(text_font, "t");	
    gStyle->SetTitleFont(text_font,"XYZ");
	gStyle->SetTitleFont(text_font,"t");  

	gStyle->SetTitleFontSize(text_size);
    gStyle->SetTextSize(text_size);	
    gStyle->SetLabelSize(text_size, "XYZ"); 
    gStyle->SetLabelSize(text_size, "t"); 	 
    gStyle->SetTitleSize(text_size, "XYZ");
    gStyle->SetTitleSize(text_size, "t");	 

    gStyle->SetTitleOffset(1., "X");
    gStyle->SetTitleOffset(1.2, "Y");  

    gStyle->SetNdivisions(ndivs, "XYZ");

	gROOT->ForceStyle();

}

//---------------------------//

// Returns bin index for value x
// Bins: [edges[i], edges[i+1])
// Last bin includes upper edge
template <size_t N>
int FindBin(const double (&edges)[N], double x)
{
    constexpr int nbins = N - 1;

    // Underflow
    if (x < edges[0])
        return -1;

    // Overflow
    if (x > edges[nbins])
        return -1;

    for (int i = 0; i < nbins; i++) {

        // Normal bins
        if (x >= edges[i] && x < edges[i+1])
            return i;
    }

    // Include right edge in last bin
    if (x == edges[nbins])
        return nbins - 1;

    return -1;
}

//---------------------------//

std::vector<TString> get_th1d_names(TFile* f) {

    std::vector<TString> histNames;

    TIter next(f->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

        TClass* cl = gROOT->GetClass(key->GetClassName());
        if (!cl) continue;

        if (cl->InheritsFrom(TH1D::Class())) {
            histNames.emplace_back(key->GetName());
        }
    }

    return histNames;
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

void calc_chi2(TH1D* h_model,
               TH1D* h_data,
               TH2D* cov,
               double &chi,
               int &ndof,
               double &pval)
{
    int NBins = cov->GetNbinsX();

    // ----------------------------
    // Build covariance matrix
    // ----------------------------
    TMatrixD V(NBins, NBins);

    for (int i = 0; i < NBins; ++i)
        for (int j = 0; j < NBins; ++j)
            V(i,j) = cov->GetBinContent(i+1, j+1);

    // Optional tiny regularization (very useful in practice)
    for (int i = 0; i < NBins; ++i)
        V(i,i) += 1e-14;

    // ----------------------------
    // Build difference vector Δ
    // ----------------------------
    TVectorD delta(NBins);

    for (int i = 0; i < NBins; ++i)
        delta[i] = h_data->GetBinContent(i+1)
                 - h_model->GetBinContent(i+1);

    // ----------------------------
    // Solve V x = Δ using Cholesky
    // ----------------------------
    TDecompChol decomp(V);

    if (!decomp.Decompose()) {
        std::cout << "ERROR: Covariance matrix not positive definite!" << std::endl;
        chi = -1;
        ndof = 0;
        pval = 0;
        return;
    }

    Bool_t ok;
    TVectorD x = decomp.Solve(delta, ok);

    if (!ok) {
        std::cout << "ERROR: Failed to solve linear system!" << std::endl;
        chi = -1;
        ndof = 0;
        pval = 0;
        return;
    }

    // ----------------------------
    // χ² = Δᵀ x
    // ----------------------------
    chi = delta * x;

    ndof = NBins;
    pval = TMath::Prob(chi, ndof);
}

//----------------------------------------//	

int ReturnIndex(double value, std::vector<double> vec) {

	int length = vec.size();
	int index = -1;

	for (int i = 0; i < length-1; i ++) {

		if (i == 0 && value == vec.at(0)) { return 0; }
		if (value > vec.at(i) && value <= vec.at(i+1)) { return i; }

	}	

	return index;

}

//----------------------------------------//

int ReturnIndexIn3DList(std::vector< std::vector< std::vector<double> > > BinEdgeVector, int FirstSliceIndex, int SecondSliceIndex, double ValueInSlice) { 

	int BinIndex = 1; // TH1D bin index, thus starting from 1
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int VectorColumnSize = BinEdgeVector.at(irow).size();

		for (int icolumn = 0; icolumn < VectorColumnSize; icolumn++){

			if (irow != FirstSliceIndex || icolumn != SecondSliceIndex) {

				BinIndex += BinEdgeVector.at(irow).at(icolumn).size()-1;

			} else {

				int LocalBins = BinEdgeVector.at(irow).at(icolumn).size();
				BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector.at(irow).at(icolumn));
				return BinIndex;

			}

		}	

	}

	return BinIndex+1; // Offset to account for bin number vs array index

}

//----------------------------------------//

std::vector<double> Return3DBinIndices(std::vector< std::vector< std::vector<double> > > BinEdgeVector) { 

	int BinCounter = 0;
	int VectorRowSize = BinEdgeVector.size();
	std::vector<double> BinIndices;

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int ielement = 0; ielement < NElements; ielement++) {

			int NElementsColumn = BinEdgeVector.at(irow).at(ielement).size();	

			for (int icolumn = 0; icolumn < NElementsColumn-1; icolumn++) {

				// Lower bin edges in the form of indices
				// + 0.5 so that the bins are centered at an integer (e.g. Bin 1, 2, 3 et al)
				BinIndices.push_back(BinCounter+0.5);
				BinCounter++;

			}	

		}

	}
	// Upper bin edge
	BinIndices.push_back(BinCounter+0.5);
	return BinIndices;

}	

//----------------------------------------//

int Return3DNBins(std::vector< std::vector< std::vector<double> > > BinEdgeVector) { 

	int NBins = 0;
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int icolumn = 0; icolumn < NElements; icolumn++) {

			int NElementsColumn = BinEdgeVector.at(irow).at(icolumn).size();

			// Number of bins for each subvector
			NBins += NElementsColumn-1;

		}

	}

	return NBins;

}

//----------------------------------------//

int ReturnIndexIn2DList(std::vector< std::vector<double> > BinEdgeVector, int SliceIndex, double ValueInSlice) { 

	int BinIndex = 1; // TH1D bin index, thus starting from 1
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		if (irow != SliceIndex) {

			BinIndex += BinEdgeVector.at(irow).size()-1;

		} else {

			int LocalBins = BinEdgeVector.at(irow).size();
			BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector.at(irow));
			return BinIndex;

		}


	}

	return BinIndex+1; // Offset to account for bin number vs array index

}

//----------------------------------------//

std::vector<double> Return2DBinIndices(std::vector< std::vector<double> > BinEdgeVector) { 

	int BinCounter = 0;
	int VectorRowSize = BinEdgeVector.size();
	std::vector<double> BinIndices;

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int ielement = 0; ielement < NElements-1; ielement++) {

			// Lower bin edges in the form of indices
			// + 0.5 so that the bins are centered at an integer (e.g. Bin 1, 2, 3 et al)
			BinIndices.push_back(BinCounter+0.5);
			BinCounter++;

		}

	}

	// Upper bin edge
	BinIndices.push_back(BinCounter+0.5);
	return BinIndices;

}	

//----------------------------------------//

int Return2DNBins(std::vector< std::vector<double> > BinEdgeVector) { 

	int NBins = 0;
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		// Number of bins for each subvector
		NBins += NElements-1;

	}

	return NBins;

}	

//---------------------------//

//----------------------------------------//
// 4D versions
//----------------------------------------//

int ReturnIndexIn4DList(const std::vector<std::vector<std::vector<std::vector<double>>>>& BinEdgeVector,
                        size_t idx0, size_t idx1, size_t idx2, double ValueInSlice) 
{
    int BinIndex = 1; // TH1D starts at 1

    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++) {
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++) {
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++) {
                if (i0 != idx0 || i1 != idx1 || i2 != idx2) {
                    BinIndex += static_cast<int>(BinEdgeVector[i0][i1][i2].size()) - 1;
                } else {
                    BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector[i0][i1][i2]);
                    return BinIndex;
                }
            }
        }
    }
    return BinIndex;
}

std::vector<double> Return4DBinIndices(const std::vector<std::vector<std::vector<std::vector<double>>>>& BinEdgeVector) {
    int BinCounter = 0;
    std::vector<double> BinIndices;

    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++) {
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++) {
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++) {
                for (size_t i3 = 0; i3 < BinEdgeVector[i0][i1][i2].size() - 1; i3++) {
                    BinIndices.push_back(BinCounter + 0.5);
                    BinCounter++;
                }
            }
        }
    }

    BinIndices.push_back(BinCounter + 0.5); // upper edge
    return BinIndices;
}

int Return4DNBins(const std::vector<std::vector<std::vector<std::vector<double>>>>& BinEdgeVector) {
    int NBins = 0;
    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++)
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++)
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++)
                NBins += static_cast<int>(BinEdgeVector[i0][i1][i2].size()) - 1;
    return NBins;
}

//----------------------------------------//
// 5D versions
//----------------------------------------//

int ReturnIndexIn5DList(const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& BinEdgeVector,
                        size_t idx0, size_t idx1, size_t idx2, size_t idx3, double ValueInSlice) 
{
    int BinIndex = 1; // TH1D starts at 1

    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++) {
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++) {
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++) {
                for (size_t i3 = 0; i3 < BinEdgeVector[i0][i1][i2].size(); i3++) {
                    if (i0 != idx0 || i1 != idx1 || i2 != idx2 || i3 != idx3) {
                        BinIndex += static_cast<int>(BinEdgeVector[i0][i1][i2][i3].size()) - 1;
                    } else {
                        BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector[i0][i1][i2][i3]);
                        return BinIndex;
                    }
                }
            }
        }
    }
    return BinIndex;
}

std::vector<double> Return5DBinIndices(const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& BinEdgeVector) {
    int BinCounter = 0;
    std::vector<double> BinIndices;

    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++)
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++)
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++)
                for (size_t i3 = 0; i3 < BinEdgeVector[i0][i1][i2].size(); i3++)
                    for (size_t i4 = 0; i4 < BinEdgeVector[i0][i1][i2][i3].size() - 1; i4++) {
                        BinIndices.push_back(BinCounter + 0.5);
                        BinCounter++;
                    }

    BinIndices.push_back(BinCounter + 0.5); // upper edge
    return BinIndices;
}

int Return5DNBins(const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& BinEdgeVector) {
    int NBins = 0;
    for (size_t i0 = 0; i0 < BinEdgeVector.size(); i0++)
        for (size_t i1 = 0; i1 < BinEdgeVector[i0].size(); i1++)
            for (size_t i2 = 0; i2 < BinEdgeVector[i0][i1].size(); i2++)
                for (size_t i3 = 0; i3 < BinEdgeVector[i0][i1][i2].size(); i3++)
                    NBins += static_cast<int>(BinEdgeVector[i0][i1][i2][i3].size()) - 1;
    return NBins;
}

//---------------------------//