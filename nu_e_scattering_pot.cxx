#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>

#include <vector>
#include <cmath>

#include "helper_functions.h"

void nu_e_scattering_pot() {

    force_plot_style();

    const double pot_ref = 1e21;
    const double events_ref = 529.0;

    const int npoints = 100;
    const double pot_max = 1.6e21;

    double efficiency = 0.4;

    std::vector<double> pot(npoints);
    std::vector<double> events(npoints);
    std::vector<double> err_x(npoints,0.0);
    std::vector<double> err_y(npoints);

    std::vector<double> events_eff(npoints);
    std::vector<double> err_eff(npoints);

    std::vector<double> frac_unc(npoints);
    std::vector<double> frac_unc_eff(npoints);

    for(int i=0;i<npoints;i++){

        pot[i] = pot_max * (double)i/(npoints-1);

        events[i] = events_ref*(pot[i]/pot_ref);
        err_y[i]  = sqrt(events[i]);

        events_eff[i] = efficiency*events[i];
        err_eff[i]    = efficiency*err_y[i];

        if(events[i]>0) frac_unc[i] = 1.0/sqrt(events[i]) * 100.;
        else frac_unc[i] = 99999;

        if(events_eff[i]>0) frac_unc_eff[i] = 1.0/sqrt(events_eff[i])*100.;
        else frac_unc_eff[i] = 99999;

    }

    double ymax = events.back()*1.15;

    TCanvas* c = new TCanvas("c","Nu-e Events vs POT",800,800);

    TPad* p1 = new TPad("p1","",0,0.3,1,1);
    TPad* p2 = new TPad("p2","",0,0,1,0.3);

    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(0.11);

    p2->SetTopMargin(0.03);
    p2->SetBottomMargin(0.35);
    p2->SetLeftMargin(0.11);

    p1->Draw();
    p2->Draw();

    //----------------------------------
    // TOP PAD
    //----------------------------------

    p1->cd();

    TGraphErrors* band = new TGraphErrors(
        npoints,&pot[0],&events[0],&err_x[0],&err_y[0]);

    band->SetFillColorAlpha(kAzure+1,0.35);
    band->SetLineColor(kAzure+2);
    band->SetLineWidth(2);

    band->SetTitle(";POT;#nu-e events");

    band->SetMaximum(ymax);
    band->SetMinimum(0);

    band->GetXaxis()->SetRangeUser(0,pot_max);
    band->GetXaxis()->SetMaxDigits(2);
    band->GetXaxis()->SetNdivisions(9);
    band->GetXaxis()->SetTitleSize(0.);
    band->GetXaxis()->SetLabelSize(0.);   

    band->GetYaxis()->CenterTitle();

    TGraphErrors* band_eff = new TGraphErrors(
        npoints,&pot[0],&events_eff[0],&err_x[0],&err_eff[0]);

    band_eff->SetFillColorAlpha(kGreen+2,0.35);
    band_eff->SetLineColor(kGreen+2);
    band_eff->SetLineWidth(2);

    band->Draw("A3L");
    band_eff->Draw("3L SAME");

    //----------------------------------
    // POT reference lines
    //----------------------------------

    double lines[] = {1e20,3e20,5.5e20,8e20};

    const char* labels[] = {
        "1#times10^{20}",
        "3#times10^{20}",
        "5.5#times10^{20}",
        "8#times10^{20}"
    };

    TLatex latex;
    latex.SetTextAlign(23);
    latex.SetTextSize(0.04);

    for(int i=0;i<4;i++){

        TLine* l = new TLine(lines[i],0,lines[i],ymax);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        //l->Draw();

        //latex.DrawLatex(lines[i],ymax*1.05,labels[i]);
    }

    //----------------------------------
    // Legend
    //----------------------------------

    TLegend* leg = new TLegend(0.55,0.75,0.85,0.88);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(132);

    leg->AddEntry(band,"100% eff","lf");
    leg->AddEntry(band_eff,"40% eff","lf");

    leg->Draw();

    //----------------------------------
    // BOTTOM PAD
    //----------------------------------

    p2->cd();

    TGraph* g_unc = new TGraph(npoints,&pot[0],&frac_unc[0]);
    TGraph* g_unc_eff = new TGraph(npoints,&pot[0],&frac_unc_eff[0]);

    g_unc->SetLineColor(kAzure+2);
    g_unc->SetLineWidth(3);

    g_unc_eff->SetLineColor(kGreen+2);
    g_unc_eff->SetLineWidth(3);

    g_unc->SetTitle(";Protons On Target (POT);uncertainty [%]");

    g_unc->SetMinimum(0);
    g_unc->SetMaximum(25);

    g_unc->GetXaxis()->SetMaxDigits(2);

    g_unc->GetXaxis()->CenterTitle();
    g_unc->GetXaxis()->SetTitleSize(0.11);
    g_unc->GetXaxis()->SetLabelSize(0.11);

    g_unc->GetXaxis()->SetRangeUser(0,pot_max);    
    g_unc->GetYaxis()->SetNdivisions(4);    
    g_unc->GetYaxis()->CenterTitle();    
    g_unc->GetYaxis()->SetTitleSize(0.11);
    g_unc->GetYaxis()->SetLabelSize(0.11);
    g_unc->GetYaxis()->SetTitleOffset(0.5);

    g_unc->Draw("AL");
    g_unc_eff->Draw("L SAME");
    //----------------------------------
    // POT reference lines
    //----------------------------------

    for(int i=0;i<4;i++){

        TLine* l = new TLine(lines[i],0,lines[i],ymax);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        //l->Draw();
        
    }


    //----------------------------------
    // Horizontal uncertainty guides
    //----------------------------------

    double guides[] = {5,10,20};

    for(int i=0;i<3;i++){

        TLine* l = new TLine(0,guides[i],pot_max,guides[i]);
        l->SetLineStyle(2);
        l->SetLineColor(kGray+1);
        l->Draw();
    }

    c->SaveAs("nu_e_scattering_pot.pdf");

    //----------------------------------
    // Extract specific POT points (40% eff)
    //----------------------------------

    auto print_point = [&](double pot_val){

        double N_eff = efficiency * events_ref * (pot_val / pot_ref);
        double stat_unc = (N_eff > 0) ? 100.0 / sqrt(N_eff) : 99999.0;

        printf("POT = %.2e --> N_eff = %.2f, stat unc = %.2f%%\n",
            pot_val, N_eff, stat_unc);
    };

    // Points of interest
    print_point(7e20);
    print_point(10.5e20);
    print_point(16e20);    
}