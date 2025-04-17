#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TChain.h"

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

int hist() {

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

string tag = "Run387440"; 


/* read in all files in the input folder */
string filename = "HIRawPrime_387440_Test.root"; 
TFile* inFile = TFile::Open(filename.c_str()); 
cout << "Opened file " << filename << endl;
TTree* zdcTree = (TTree*)inFile->Get("zdcanalyzer/zdcrechit");


/* read in the unpacked information */
TTreeReader zdcReader(zdcTree);
TTreeReaderValue<float> ZDCp_Energy(zdcReader, "ZDCp_Energy");
TTreeReaderValue<float> ZDCm_Energy(zdcReader, "ZDCm_Energy");


/* create histograms */

TH1D* hZDCp_En = new TH1D("hZDCp_En", "ZDC Energy Plus", 200, 0, 20000);
TH1D* hZDCm_En = new TH1D("hZDCm_En", "ZDC Energy Minus", 400, 0, 20000);

Long64_t totalEvents = zdcReader.GetEntries(true);
for (Long64_t i = 0; i < totalEvents; i++) {
zdcReader.Next();

if (i % 10000 == 0) {
cout << "Event " << i << endl;
}
hZDCp_En->Fill(*ZDCp_Energy);
hZDCm_En->Fill(*ZDCm_Energy);
} // end loop over  the  number of events

std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
std::cout << "Finished processing " << totalEvents << " events." << std::endl;
std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;


TFile* outfile = new TFile(Form("outfile_%s.root", tag.c_str()), "RECREATE");
outfile->cd();
hZDCp_En->Write();
hZDCm_En->Write();
outfile->Close();


TLatex* cms = new TLatex(0.10,0.92,"#bf{CMS} Run 387440");
cms->SetNDC();
cms->SetTextSize(0.05);
cms->SetTextFont(42);

// ------------------------------------------
/* ZDC Plus */
// ------------------------------------------

TLegend* leg3 = new TLegend(0.65,0.7,0.8,0.85);
leg3->SetBorderSize(0);
leg3->SetFillStyle(0);
leg3->SetTextSize(0.045);

TCanvas* c3 = new TCanvas("c3","c3",800,600);
c3->cd();
c3->SetTickx(1);
c3->SetTicky(1);
c3->SetLogy(); 
c3->SetTopMargin(0.09);
c3->SetBottomMargin(0.11);
c3->SetLeftMargin(0.09);
c3->SetRightMargin(0.05);

// zdc plus
hZDCp_En->SetLineColor(kRed);
hZDCp_En->SetLineWidth(2);
hZDCp_En->GetXaxis()->SetTitleSize(0.05);
hZDCp_En->SetMarkerStyle(20);
hZDCp_En->GetXaxis()->SetTitle("ZDC Energy(GeV)");
hZDCm_En->SetLineColor(kBlue);
hZDCm_En->SetLineWidth(2);
hZDCm_En->GetXaxis()->SetTitleSize(0.05);
hZDCm_En->SetMarkerStyle(20);
hZDCm_En->GetXaxis()->SetTitle("ZDC Energy(GeV)");
hZDCp_En->Draw();
hZDCm_En->Draw("same");
leg3->AddEntry(hZDCp_En,"ZDC Plus","l");
leg3->AddEntry(hZDCm_En,"ZDC Minus","l");
leg3->Draw("same");
cms->Draw("same");
c3->SaveAs(Form("ZDC_En%s.pdf", tag.c_str()));

return 1;
}

int main(int argc, char *argv[]) {
if (argc == 2)
return hist();
else {
cout << "ERROR" << endl;
return -1;
}
}