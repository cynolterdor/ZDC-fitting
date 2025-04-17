#include "TFile.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"

#include "TCanvas.h"
#include "TLegend.h"

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

//function to find the bin with the minimum content in a given range
double min_finder(int low, int high, TH1D* hist) {
    int bin_min = hist->FindBin(low);
    int bin_max = hist->FindBin(high);
    double min_content = 1e9;
    int min_bin = 0;
    for (int i = bin_min; i <= bin_max; i++) {
        double bin_content = hist->GetBinContent(i);
        if (bin_content < min_content) {
            min_content = bin_content;
            min_bin = i;
        }
    }
    double x_fit = hist->GetXaxis()->GetBinUpEdge(min_bin);
    return x_fit;
}

int fitter() {

    //hide statistics box and title
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    //run no
    string tag = "Run387440"; 

    string filename = "outfile_Run387440.root"; 
    TFile* inFile = TFile::Open(filename.c_str()); 
    cout << "Opened file " << filename << endl;
    TH1D* hZDCp_En = (TH1D*)inFile->Get("hZDCp_En");
    TH1D* hZDCm_En = (TH1D*)inFile->Get("hZDCm_En");




    // ------------------------------------------
    /* ZDC Plus */
    // ------------------------------------------
    TCanvas *c1 = new TCanvas("c1","Fit",800,600);
    c1->cd();
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLogy(); 
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.11);
    c1->SetLeftMargin(0.09);
    c1->SetRightMargin(0.05);
    
    //fitting exponential to background
    int x_low = 1;
    int x_high = 1000;
    double x_fit1p = min_finder(x_low, x_high, hZDCp_En);
    cout << "Fitting exponential function to ZDC Plus data from 0 to " << x_fit1p << endl;
    TF1 *f1 = new TF1("f1","exp([0]+[1]*x)/(x-[2])", 0, 3000);
    hZDCp_En->Fit("f1","R","same",0,x_fit1p);

    //fitting gaussian to first peak
    x_low = 5000;
    x_high = 6000;
    double x_fit2p = min_finder(x_low, x_high, hZDCp_En);
    cout << "Fitting gaussian function to ZDC Plus data from " << x_fit1p << " to " << x_fit2p << endl;
    TF1 *f2 = new TF1("f2","gaus", 0, 8000);
    hZDCp_En->Fit("f2","R","same",x_fit1p,x_fit2p);

    //fitting gaussian to second peak
    x_low = 8000;
    x_high = 9500;
    double x_fit3p = min_finder(x_low, x_high, hZDCp_En);
    cout << "Fitting gaussian function to ZDC Plus data from " << x_fit2p << " to " << x_fit3p << endl;
    TF1 *f3 = new TF1("f3","gaus", 0, 16000);
    hZDCp_En->Fit("f3","R","same",x_fit2p,x_fit3p);

    //all fits together
    TF1 *f_all_p = new TF1("f_all_p", "exp([0] + [1]*x)/(x-[2]) + [3]*exp(-0.5*((x-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2)", 0, 10000);
    f_all_p->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f2->GetParameter(0), f2->GetParameter(1), f2->GetParameter(2), f3->GetParameter(0), f3->GetParameter(1), f3->GetParameter(2));
    hZDCp_En->Fit("f1","R","same",0,x_fit1p);
    cout << "Fitting function to ZDC Plus data from 0 to 10000" << endl;
    hZDCp_En->Fit("f_all_p", "R"  , "same", 0 , 10000);

    f1->SetParameters(f_all_p->GetParameter(0), f_all_p->GetParameter(1), f_all_p->GetParameter(2));
    f2->SetParameters(f_all_p->GetParameter(3), f_all_p->GetParameter(4), f_all_p->GetParameter(5));
    f3->SetParameters(f_all_p->GetParameter(6), f_all_p->GetParameter(7), f_all_p->GetParameter(8));
    
    //efficiency function
    TF1 *f1_erf = new TF1("f1_erf", "1 - [0]*(1 + TMath::Erf((x-[1])/sqrt(2)/[2]))", 0, 6000);
    f1_erf->SetParameters(0.5, f2->GetParameter(1), f2->GetParameter(2));


     

    //Make plot pretty
    TLatex* cms = new TLatex(0.10,0.92,"#bf{CMS} Run 387440");
    cms->SetNDC();
    cms->SetTextSize(0.05);
    cms->SetTextFont(42);

    TLatex* title_m = new TLatex(0.7,0.92,"ZDC Minus");
    title_m->SetNDC();
    title_m->SetTextSize(0.05);
    title_m->SetTextFont(42);

    TLatex* title_p = new TLatex(0.7,0.92,"ZDC Plus");
    title_p->SetNDC();
    title_p->SetTextSize(0.05);
    title_p->SetTextFont(42);

    TLegend* leg1 = new TLegend(0.5,0.7,0.8,0.85);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.045);
    TLegend* leg2 = new TLegend(0.65,0.7,0.8,0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.045);

    hZDCp_En->SetLineColor(kBlue);
    hZDCp_En->SetLineWidth(2);
    hZDCp_En->GetXaxis()->SetTitleSize(0.05);
    hZDCp_En->SetMarkerStyle(20);
    hZDCp_En->GetXaxis()->SetTitle("ZDC Energy(GeV)");

    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);


    //Draw the histogram and the fits
    hZDCp_En->GetXaxis()->SetRangeUser(0, 10000);
    hZDCp_En->Draw();
    f2->SetRange(0, 10000);
    f2->SetLineColor(kGray);
    f2->SetLineWidth(2);
    f2->SetLineStyle(9);
    f2->Draw("same");
;
    f3->SetRange(0, 10000);
    f3->SetLineColor(kGray);
    f3->SetLineWidth(2);
    f3->SetLineStyle(9);
    f3->Draw("same");

    f_all_p->Draw("same");

    hZDCp_En->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    leg1->AddEntry(hZDCp_En,"ZDC Plus","l");
    leg1->AddEntry(f1,"Fit","l");
    leg1->AddEntry(f2,"Extrapolated Gaussians","l");
    leg1->Draw("same");
    cms->Draw("same");
    title_p->Draw("same");
    c1->SaveAs(Form("Fitted_%s_p.pdf", tag.c_str()));

    TCanvas *c2 = new TCanvas("c2","Cut_eff",800,600);
    c2->cd();
    c2->SetTickx(1);
    c2->SetTicky(1);
    c2->SetTopMargin(0.09);
    c2->SetBottomMargin(0.11);
    c2->SetLeftMargin(0.09);
    c2->SetRightMargin(0.05);

    f1_erf->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    f1_erf->GetYaxis()->SetTitle("Efficiency");
    f1_erf->Draw();
    cms->Draw("same");
    title_p->Draw("same");
    c2->SaveAs(Form("Cut_eff_%s_p.pdf", tag.c_str()));

    TCanvas *c3 = new TCanvas("c3","Pull_plus",800,600);
    c3->cd();
    c3->SetTickx(1);
    c3->SetTicky(1);
    c3->SetTopMargin(0.09);
    c3->SetBottomMargin(0.11);
    c3->SetLeftMargin(0.09);
    c3->SetRightMargin(0.05);

    hZDCp_En->Sumw2();

    TH1F *hPull_p = (TH1F*)hZDCp_En->Clone("hPull_p");
    TH1F *hPull_p1 = (TH1F*)hZDCp_En->Clone("hPull_p1");
    TH1F *hPull_p2 = (TH1F*)hZDCp_En->Clone("hPull_p2");
    hPull_p->Reset();
    hPull_p1->Reset();
    hPull_p2->Reset();

    /*for (int i = 1; i <= hZDCp_En->GetNbinsX(); i++) {
        double data_p = hZDCp_En->GetBinContent(i);
        //double fit_p = f1_3->Eval(hZDCp_En->GetBinCenter(i));
        //double fit_p1 = f1_1->Eval(hZDCp_En->GetBinCenter(i));
        //double fit_p2 = f1->Eval(hZDCp_En->GetBinCenter(i));
        double error_p = hZDCp_En->GetBinError(i);
        
        if (error_p > 0) {
            hPull_p->SetBinContent(i, (fit_p - data_p) / error_p);
            hPull_p1->SetBinContent(i, (fit_p1 - data_p) / error_p);
            hPull_p2->SetBinContent(i, (fit_p2 - data_p) / error_p);
        }
    }*/
    

    /*hPull_p->GetXaxis()->SetRangeUser(0, x_fit1p);
    //hPull_p->GetYaxis()->SetRangeUser(-0.5, 0.5);
    hPull_p->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    hPull_p->GetYaxis()->SetTitle("Pull");
    hPull_p->SetLineColor(kOrange);
    hPull_p->Draw("HIST");
    //hPull_p->GetXaxis()->SetRangeUser(x_fit1p_1, x_fit1p_2);
    hPull_p1->SetLineColor(kGreen);
    hPull_p1->Draw("HIST SAME");
    //hPull_p->GetXaxis()->SetRangeUser(x_fit1p_2, x_fit1p);
    hPull_p2->SetLineColor(kRed);
    hPull_p2->Draw("HIST SAME");
    cms->Draw("same");
    title_p->Draw("same");
    c3->SaveAs(Form("Pull_plus.pdf"));*/



    // ------------------------------------------
    /* ZDC Minus */
    // ------------------------------------------
    TCanvas *c4 = new TCanvas("c4","fit_minus",800,600);
    c4->cd();
    c4->SetTickx(1);
    c4->SetTicky(1);
    c4->SetLogy();
    c4->SetTopMargin(0.09);
    c4->SetBottomMargin(0.11);
    c4->SetLeftMargin(0.09);
    c4->SetRightMargin(0.05);
    

    //fitting exponential to background  
    x_low = 0;
    x_high = 1000;
    double x_fit1m = min_finder(x_low, x_high, hZDCm_En);
    cout << "Fitting exponential function to ZDC Minus data from 0 to " << x_fit1m << endl;
    TF1 *f4 = new TF1("f4","exp([0] + [1]*x)/(x-[2])**4 + [3]/(x-[4])**3", 0, 3000);
    hZDCm_En->Fit("f4","R","same",0,x_fit1m);

    //fit gaussian to first peak
    x_low = 3000;
    x_high = 4000;
    double x_fit2m = min_finder(x_low, x_high, hZDCm_En);
    cout << "Fitting gaussian function to ZDC Minus data from " << x_fit1m << " to " << x_fit2m << endl;
    TF1 *f5 = new TF1("f5","gaus", 0, 8000);
    hZDCm_En->Fit("f5","R","same",1000,3000);


    //fitting gaussian to second peak
    x_low = 5000;
    x_high = 6000;
    double x_fit3m = min_finder(x_low, x_high, hZDCm_En);
    cout << "Fitting gaussian function to ZDC Minus data from " << x_fit2m << " to " << x_fit3m << endl;
    TF1 *f6 = new TF1("f6","gaus", 0, 16000);
    hZDCm_En->Fit("f6","R","same",x_fit2m,x_fit3m);

    //all fits together
    TF1 *f_all_m = new TF1("f_all_m", "exp([0] + [1]*x)/(x-[2])**4 + [3]/(x-[4])**3 +[5]*exp(-0.5*((x-[6])/[7])**2) + [8]*exp(-0.5*((x-[9])/[10])**2)", 0, 10000);
    f_all_m->SetParameters(f4->GetParameter(0), f4->GetParameter(1),f4->GetParameter(2),f4->GetParameter(3), f4->GetParameter(4), f5->GetParameter(0), f5->GetParameter(1), f5->GetParameter(2), f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2));
    cout << "Fitting function to ZDC Minus data from 0 to 10000" << endl;
    hZDCm_En->Fit("f_all_m","R","same",0,7000);

    f4->SetParameters(f_all_m->GetParameter(0), f_all_m->GetParameter(1), f_all_m->GetParameter(2)), f_all_m->GetParameter(3), f_all_m->GetParameter(4);
    f5->SetParameters(f_all_m->GetParameter(5), f_all_m->GetParameter(6), f_all_m->GetParameter(7));
    f6->SetParameters(f_all_m->GetParameter(8), f_all_m->GetParameter(9), f_all_m->GetParameter(10));


    //efficiency function
    TF1 *f2_erf = new TF1("f2_erf", "1 - [0]*(1 + TMath::Erf((x-[1])/sqrt(2)/[2]))", 0, 4000);
    f2_erf->SetParameters(0.5, f5->GetParameter(1), f5->GetParameter(2));


    //Draw the histogram and the fits
    hZDCm_En->GetXaxis()->SetRangeUser(0, 7000);
    hZDCm_En->Draw();

    f4->SetLineColor(kGray);
    f4->SetLineStyle(9);
    f4->Draw("same");
    f5->SetLineColor(kGray);
    f5->SetLineStyle(9);
    f5->Draw("same");
    f6->SetLineColor(kGray);
    f6->SetLineStyle(9);
    f6->Draw("same");
    f_all_m->SetLineColor(kRed);
    f_all_m->Draw("same");


    hZDCm_En->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    leg2->AddEntry(hZDCm_En,"ZDC Minus","l");
    leg2->AddEntry(f_all_m,"Fit","l");
    leg2->AddEntry(f5,"Extrapolated fits","l");
    leg2->Draw("same");
    cms->Draw("same");
    title_m->Draw("same");
    c4->SaveAs(Form("Fitted_%s_m.pdf", tag.c_str()));

    TCanvas *c5 = new TCanvas("c5","Cut_eff",800,600);
    c5->cd();
    c5->SetTickx(1);
    c5->SetTicky(1);

    f2_erf->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    f2_erf->GetYaxis()->SetTitle("Efficiency");
    f2_erf->Draw();
    cms->Draw("same");
    title_m->Draw("same");
    c5->SaveAs(Form("Cut_eff_%s_m.pdf", tag.c_str()));
    
    /*TCanvas *c6 = new TCanvas("c6", "Pull", 800, 600);
    c6->cd();
    c6->SetTickx(1);
    c6->SetTicky(1);
    c6->SetTopMargin(0.09);
    c6->SetBottomMargin(0.11);
    c6->SetLeftMargin(0.09);
    c6->SetRightMargin(0.05);

    hZDCm_En->Sumw2();

    TH1F *hPull = (TH1F*)hZDCm_En->Clone("hPull");
    hPull->Reset();

    for (int i = 1; i <= hZDCm_En->GetNbinsX(); i++) {
        double data = hZDCm_En->GetBinContent(i);
        double fit = f4_2->Eval(hZDCm_En->GetBinCenter(i));
        double error = hZDCm_En->GetBinError(i);
        
        if (error > 0) hPull->SetBinContent(i, (fit - data) / error);
    }

    hPull->GetXaxis()->SetRangeUser(0, x_fit1m);
    hPull->GetXaxis()->SetTitle("ZDC Energy(GeV)");
    hPull->GetYaxis()->SetTitle("Pull");
    hPull->Draw();
    cms->Draw("same");
    title_m->Draw("same");
    c6->SaveAs(Form("Pull_minus.pdf"));*/
    
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc == 2)
    return fitter();
    else {
    cout << "ERROR" << endl;
    return -1;
    }
}