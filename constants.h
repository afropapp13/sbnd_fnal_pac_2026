#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <array>
#include <cstddef> // std::size

using namespace std;

namespace constants {

	//----------------------------------------//

	static const double A = 40.;

	static const double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}    

    const int text_font = 132;
    const double text_size = 0.06;
    const int ndivs = 6; 

    double ProtonMass_GeV = 0.9382720813;

	//----------------------------------------//

    static constexpr double ArrayNBinsDeltaPT[] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9};
    static constexpr double ArrayNBinsDeltaAlphaT[] = {0.,22.,44.,66.,88.,110.,145.,180.};
    static constexpr double ArrayNBinsMuonMomentum[] = {0.1,0.2,0.3,0.4,0.5,0.64,0.77,0.9,1.,1.1,1.2};
    static constexpr double ArrayNBinsMuonCosTheta[] = {-1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

    // Number of bins = number of edges - 1
    static constexpr std::size_t NBinsDeltaPT = std::size(ArrayNBinsDeltaPT) - 1;
    static constexpr std::size_t NBinsDeltaAlphaT = std::size(ArrayNBinsDeltaAlphaT) - 1;
    static constexpr std::size_t NBinsMuonMomentum = std::size(ArrayNBinsMuonMomentum) - 1;
    static constexpr std::size_t NBinsMuonCosTheta = std::size(ArrayNBinsMuonCosTheta) - 1;

	//----------------------------------------//
	
	// And moving on towards a 2D analysis
		
	static const int TwoDNBinsDeltaPT = 3; std::vector<double> TwoDArrayNBinsDeltaPT{0.0,0.2,0.4,1.0};
	static const int TwoDNBinsDeltaAlphaT = 4; std::vector<double> TwoDArrayNBinsDeltaAlphaT{0.0,45.0,90.0,135.0,180.0};
	static const int TwoDNBinsMuonCosTheta = 4; std::vector<double> TwoDArrayNBinsMuonCosTheta{-1.,0.0,0.5,0.75,1.0};
	static const int TwoDNBinsMuonMomentum = 3; std::vector<double> TwoDArrayNBinsMuonMomentum{0.1,0.4,0.6,1.2};
	
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.}, // DeltaPT < 0.2 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.2 < DeltaPT < 0.4 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}  // DeltaPT > 0.4 GeV/c
												
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices{ 
												{ 0.1,0.2,0.3,0.4,1.2},                           // -1    < cosθμ < 0
												{ 0.1,0.2,0.3,0.4,0.5,1.2},                       //  0    < cosθμ < 0.5
												{ 0.1,0.2,0.3,0.4,0.5,0.64,0.77,1.2},             //  0.5  < cosθμ < 0.75
												{ 0.1,0.2,0.3,0.4,0.5,0.64,0.77,0.9,1.,1.1,1.2}   //  0.75 < cosθμ < 1
											};
											
	//----------------------------------------//	
	
	std::vector< std::vector< std::vector<double> > > TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices{
														{ // -1 < cosθμ < 0
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,1.6},      // 0.40 < Pμ < 0.60 GeV/c
														{0.2,1.6}           // 0.60 < Pμ < 1.20 GeV/c
														},	
														
														{ // 0 < cosθμ < 0.5
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,0.8,1.6},  // 0.40 < Pμ < 0.6 GeV/c
														{0.2,1.6}   // 0.6 < Pμ < 1.20 GeV/c
														},
														
														{ // 0.5 < cosθμ < 0.75
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,0.8,1.6},  // 0.40 < Pμ < 0.6 GeV/c
														{0.2,0.8,1.,1.6}    // 0.6 < Pμ < 1.20 GeV/c
														},														
														
														{ // 0.75 < cosθμ < 1
														{0.2,0.5,1.6},      // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.8,1.6},      // 0.40 < Pμ < 0.6 GeV/c
														{0.2,0.8,1.2,1.6}   // 0.6 < Pμ < 1.20 GeV/c
														}		
																
													};
	//----------------------------------------//													

	// Scaling factor for multi dimensional analysis

	static std::map<TString,double> MultiDimScaleFactor =
	{
		{ "DeltaPTPlot", 1. },
		{ "DeltaAlphaTPlot", 1. },		
		{ "MuonMomentumPlot", 1. },
		{ "MuonCosThetaPlot", 1. },
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot", TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot", TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) },
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot", TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) },
		{ "DeltaPT_DeltaAlphaT_0_00To45_00Plot", TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot", TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) },
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot", TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot", TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) },
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot", TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot", TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot", TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) },
	
    };	

}
#endif