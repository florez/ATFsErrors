#include<iostream>
#include <fstream>
#include<vector>
#include<string>
#include<TROOT.h>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TRandom.h"

using namespace std;

void ATFPlotter_usingMeanX(){

   gStyle->SetOptStat(0); // this line just removes the stats box from the plot

  // Set your input data files
  TFile *dataOS = TFile::Open("mutau_data/2019/Data_OS.root"); // data root file for OS CR
  TFile *dataLS = TFile::Open("mutau_data/2019/Data_LS.root"); // data root file for LS CR

  // Load your histos
  dataOS->cd();
  TH1D *dataOS_LargestDiTauMass = (TH1D*)gDirectory->Get("Muon1Tau1ReconstructableMassDeltaPt");
  Double_t nbins_ditaumass = dataOS_LargestDiTauMass->GetNbinsX();  //get the number of bins for the histo. for ditau mass.
  Double_t xmin_ditaumass = dataOS_LargestDiTauMass->GetXaxis()->GetXmin();  //get the minimum for the x-axis for the histo. for ditau mass.
  Double_t xmax_ditaumass = dataOS_LargestDiTauMass->GetXaxis()->GetXmax();  //get the maximum for the x-axis for the histo. for ditau mass.

  // the lines below are just making a copy of the ditau mass histogram (we don't like to use the initial clone because ROOT creates memory issues)
  TH1D *h_dataOS = new TH1D("h_dataOS", "h_dataOS", nbins_ditaumass, xmin_ditaumass, xmax_ditaumass);
  for(int i=1; i<= (nbins_ditaumass+1); i++) {
    h_dataOS->SetBinContent(i,dataOS_LargestDiTauMass->GetBinContent(i));
    h_dataOS->SetBinError(i,dataOS_LargestDiTauMass->GetBinError(i));      
  }

  // go to the OS root file, and then cd into the appropriate folder (in this case NDiJetCombinations)
  //dataLS->cd("NDiJetCombinations");
  dataLS->cd();
  TH1D *dataLS_LargestDiTauMass = (TH1D*)gDirectory->Get("Muon1Tau1ReconstructableMassDeltaPt");

  nbins_ditaumass = dataLS_LargestDiTauMass->GetNbinsX();  //get the number of bins for the histo. for ditau mass.
  xmin_ditaumass = dataLS_LargestDiTauMass->GetXaxis()->GetXmin();  //get the minimum for the x-axis for the histo. for ditau mass.
  xmax_ditaumass = dataLS_LargestDiTauMass->GetXaxis()->GetXmax();  //get the maximum for the x-axis for the histo. for ditau mass.

  // the lines below are just making a copy of the ditau mass histogram (we don't like to use the initial clone because ROOT creates memory issues)
  TH1D *h_dataLS = new TH1D("h_dataLS", "h_dataLS", nbins_ditaumass, xmin_ditaumass, xmax_ditaumass);
  for(int i=1; i<= (nbins_ditaumass); i++) {
    h_dataLS->SetBinContent(i,dataLS_LargestDiTauMass->GetBinContent(i));
    h_dataLS->SetBinError(i,dataLS_LargestDiTauMass->GetBinError(i));
  }

  // we're going to rebin the default histograms
  const int theNbins = 4;
  Double_t xbins_mass[theNbins+1] = {0, 200, 400, 600, 1000}; // this is the binning we want for the mass-dependent ATF
 
  Double_t xbins_massMeans[theNbins]; // array which will contain the means of each mass bin
  Double_t xbins_massStdDev[theNbins]; // array which will contain the means of each mass bin

  Int_t xlowrange = 0;
  Int_t  xhighrange = 0;
  Double_t binMean = 0;
  Double_t binStdDev = 0;
  for(int i=1; i<=theNbins; i++) {

    if(i == theNbins) {
      xlowrange = h_dataOS->FindBin(xbins_mass[i-1]);
      xhighrange = nbins_ditaumass + 1;
    } else {
      xlowrange = h_dataOS->FindBin(xbins_mass[i-1]);
      xhighrange = h_dataOS->FindBin(xbins_mass[i]) - 1;
    }
    
    h_dataOS->GetXaxis()->SetRange(xlowrange,xhighrange);
    binMean = h_dataOS->GetMean();
    binStdDev = h_dataOS->GetStdDev();
    xbins_massMeans[i-1] = binMean;
    // If you have low stats in your analysis, comment the line below and un-comment the line that follows.
    // This basically sets the erros on X to zero and avoid drastic fluctuations of the fit that are not
    // correlarted with the agremeent of the points w.r.t the pol0 fit. 
    xbins_massStdDev[i-1] = binStdDev;
    //xbins_massStdDev[i-1] = 0;
    std::cout << "Bin " << i << " mean = " << xbins_massMeans[i-1] << std::endl;    
  }

  // re-bin all histograms below
  TH1 *h_dataOS_rebinned = h_dataOS->Rebin(theNbins,"h_dataOS_rebinned",xbins_mass);
  TH1 *h_dataLS_rebinned = h_dataLS->Rebin(theNbins,"h_dataLS_rebinned",xbins_mass);

  h_dataOS_rebinned->SetBinContent(theNbins, h_dataOS_rebinned->GetBinContent(theNbins) + h_dataOS_rebinned->GetBinContent(theNbins+1));
  h_dataOS_rebinned->SetBinError(theNbins, sqrt(pow(h_dataOS_rebinned->GetBinError(theNbins),2.)+pow(h_dataOS_rebinned->GetBinError(theNbins+1),2.)));
  h_dataOS_rebinned->SetBinContent(theNbins+1, 0);
  h_dataOS_rebinned->SetBinError(theNbins+1, 0);

  h_dataLS_rebinned->SetBinContent(theNbins, h_dataLS_rebinned->GetBinContent(theNbins) + h_dataLS_rebinned->GetBinContent(theNbins+1));
  h_dataLS_rebinned->SetBinError(theNbins, sqrt(pow(h_dataLS_rebinned->GetBinError(theNbins),2.)+pow(h_dataLS_rebinned->GetBinError(theNbins+1),2.)));
  h_dataLS_rebinned->SetBinContent(theNbins+1, 0);
  h_dataLS_rebinned->SetBinError(theNbins+1, 0);

  std::cout << "OS yield " << h_dataOS_rebinned->Integral(0,theNbins) << std::endl;

  h_dataOS_rebinned->Divide(h_dataLS_rebinned);
  Double_t ybins_ATF[theNbins];
  Double_t ybins_ATFerrors[theNbins];

  for(int i=1; i<=theNbins; i++) {
    ybins_ATF[i-1] = h_dataOS_rebinned->GetBinContent(i);
    ybins_ATFerrors[i-1] = h_dataOS_rebinned->GetBinError(i);
  }

   gStyle->SetOptStat(0); // this line just removes the stats box from the plot

   TCanvas *c1 = new TCanvas("c1", "c1",158,150,700,500);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(3);
   c1->SetLeftMargin(0.08);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetFrameBorderMode(0);

// ------------> top pad of the plot
   TPad *top = new TPad("top", "top",0,0.25,1,1);
   top->SetTickx();
   top->SetTicky();
   top->Draw();
   top->cd();
   top->Range(61.49425,-1.836137,167.8161,16.52523);
   top->SetFillColor(0);
   top->SetBorderMode(-1);
   top->SetBorderSize(5);
   top->SetLeftMargin(0.08);
   top->SetRightMargin(0.05);
   top->SetTopMargin(0.05);
   top->SetFrameBorderMode(0);
   top->SetFrameBorderMode(0);

  TGraphErrors *ATF = new TGraphErrors(theNbins,xbins_massMeans,ybins_ATF,xbins_massStdDev,ybins_ATFerrors);
  ATF->SetMarkerColor(kBlack);
  ATF->GetYaxis()->SetTitle("ATF");
  ATF->GetYaxis()->SetTitleSize(0.06);
  ATF->GetYaxis()->SetTitleOffset(0.65);
  ATF->SetLineColor(kBlack);
  ATF->SetMarkerStyle(20);
  ATF->SetMarkerSize(1.0);
  ATF->Draw("AEP"); // draw the observed data in the opposite-sign pass-VBF region

  // The transfer factor is mass dependent, so we're going to fit the TF data points to a 1D line. We're creating this function below
  TF1 *pol0fit = new TF1("pol0fit","[0]",xbins_massMeans[0] - xbins_massStdDev[0] - 1.,xbins_massMeans[theNbins-1] + xbins_massStdDev[theNbins-1] + 1.); // [0] and [1] are the parameters of the 1D polynomial (i.e. y-intercept and slope). NOTE: fit is only defined from above 200 GeV

  ATF->Fit("pol0fit","R"); // fit the TF data points using our fit function
  double pol0fit_param0 = pol0fit->GetParameter(0); // extract the parameters of the fit function

  // The transfer factor is mass dependent, so we're going to fit the TF data points to a 1D line. We're creating this function below
  TF1 *pol1fit = new TF1("pol1fit","[0]+[1]*x",xbins_massMeans[0] - xbins_massStdDev[0] - 1.,xbins_massMeans[theNbins-1] + xbins_massStdDev[theNbins-1] +1.); // [0] and [1] are the parameters of the 1D polynomial (i.e. y-intercept and slope). NOTE: fit is only defined from above 200 GeV
  pol1fit->SetParameters(1,0); // set dummmy values to start the fit ... using a "flat line" as a starting point

  ATF->Fit("pol1fit","R"); // fit the TF data points using our fit function
  double pol1fit_param0 = pol1fit->GetParameter(0); // extract the parameters of the fit function
  double pol1fit_param1 = pol1fit->GetParameter(1); // extract the parameters of the fit function


  for(int i=1; i<=theNbins; i++) {
    std::cout << "Systematic error in bin " << i << " = " << (abs((pol1fit_param0 + (pol1fit_param1 * xbins_massMeans[i-1])) - pol0fit_param0) / pol0fit_param0 * 100.0) << "%" << std::endl;
  }

   TF1 *PrevFitTMP1 = new TF1("PrevFitTMP","pol0",xbins_massMeans[0] - xbins_massStdDev[0],xbins_massMeans[theNbins-1] + xbins_massStdDev[theNbins-1], TF1::EAddToList::kNo);
   PrevFitTMP1->SetFillColor(19);
   PrevFitTMP1->SetFillStyle(0);
   PrevFitTMP1->SetMarkerStyle(20);
   PrevFitTMP1->SetLineColor(2);
   PrevFitTMP1->SetLineWidth(1);
   PrevFitTMP1->SetLineStyle(0);
   PrevFitTMP1->GetXaxis()->SetLabelFont(42);
   PrevFitTMP1->GetXaxis()->SetLabelSize(0.035);
   PrevFitTMP1->GetXaxis()->SetTitleSize(0.05);
   PrevFitTMP1->GetXaxis()->SetTitleOffset(1);
   PrevFitTMP1->GetXaxis()->SetTitleFont(42);
   PrevFitTMP1->GetYaxis()->SetLabelFont(42);
   PrevFitTMP1->GetYaxis()->SetLabelSize(0.035);
   PrevFitTMP1->GetYaxis()->SetTitleSize(0.035);
   PrevFitTMP1->GetYaxis()->SetTitleFont(42);
   PrevFitTMP1->SetParameter(0,pol0fit_param0); // flat fit line
   PrevFitTMP1->SetParError(0,0);
   PrevFitTMP1->Draw("same");

  top->Modified();
  c1->cd();

// ------------>ratio pad, which contains the ratio of the OS-failVBF observation to the OS-failVBF prediction using the TF
  TPad *bottom = new TPad("bottom", "",0,0,1,0.32);
  bottom->SetTickx();
  bottom->SetTicky();
  bottom->Draw();
  bottom->cd();
//  bottom->Range(61.49425,-1.281439,167.8161,2.98999);
  bottom->SetFillColor(0);
  bottom->SetBorderMode(-1);
  bottom->SetBorderSize(5);
  bottom->SetLeftMargin(0.08);
  bottom->SetRightMargin(0.05);
  bottom->SetTopMargin(0);
  bottom->SetBottomMargin(0.3);
  bottom->SetFrameBorderMode(0);
  bottom->SetFrameBorderMode(0);

  // Below we're adding the uncertainty band to the ratio pad
  Double_t Graph0_fy1002[theNbins];
  Double_t Graph0_fx1002[theNbins];
  Double_t Graph0_fex1002[theNbins];
  Double_t Graph0_fey1002[theNbins];
  for(int i=1; i<= theNbins; i++) {
    Graph0_fy1002[i-1] = 0;
    Graph0_fey1002[i-1] = (abs((pol1fit_param0 + (pol1fit_param1 * xbins_massMeans[i-1])) - pol0fit_param0) / pol0fit_param0);
    Graph0_fx1002[i-1] = (xbins_mass[i] + xbins_mass[i-1]) / 2.0;
    Graph0_fex1002[i-1] = (xbins_mass[i] - xbins_mass[i-1]) / 2.0;
  }

   TGraphErrors *gre = new TGraphErrors(theNbins,Graph0_fx1002,Graph0_fy1002,Graph0_fex1002,Graph0_fey1002);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc00cc");
   gre->SetFillColor(ci);
   gre->SetFillStyle(3004);

   gre->SetMinimum(-0.5);
   gre->SetMaximum(0.5);
   gre->GetXaxis()->SetTitle("m(#tau,#tau,#slash{E}_{T}) [GeV]");
   gre->GetXaxis()->SetLabelFont(42);
   gre->GetXaxis()->SetLabelOffset(0.007);
   gre->GetXaxis()->SetLabelSize(0.105);
   gre->GetXaxis()->SetTitle("m(#tau, #tau, #Delta p_{T}) [GeV]");
   gre->GetXaxis()->SetTitleSize(0.12);
   gre->GetXaxis()->SetTitleOffset(0.9);
   gre->GetXaxis()->SetTitleFont(42);
   gre->GetYaxis()->SetTitleOffset(0.35);
   gre->GetYaxis()->SetTitle("Syst. Error");
   gre->GetYaxis()->SetNdivisions(505);
   gre->GetYaxis()->SetLabelFont(42);
   gre->GetYaxis()->SetLabelOffset(0.007);
   gre->GetYaxis()->SetLabelSize(0.105);
   gre->GetYaxis()->SetTitleSize(0.105);
   gre->GetYaxis()->SetTitleFont(42);
//   gre->Draw("AXIS");
   gre->Draw("A2"); // draw the uncertainty band around 1

   TF1 *PrevFitTMP2 = new TF1("PrevFitTMP2","pol0",xbins_mass[0],xbins_mass[theNbins], TF1::EAddToList::kNo);
   PrevFitTMP2->SetLineColor(2);
   PrevFitTMP2->SetLineWidth(1);
   PrevFitTMP2->SetLineStyle(0);
   PrevFitTMP2->GetXaxis()->SetLabelFont(42);
   PrevFitTMP2->GetXaxis()->SetLabelSize(0.035);
   PrevFitTMP2->GetXaxis()->SetTitleSize(0.035);
   PrevFitTMP2->GetXaxis()->SetTitleOffset(1);
   PrevFitTMP2->GetXaxis()->SetTitleFont(42);
   PrevFitTMP2->GetYaxis()->SetLabelFont(42);
   PrevFitTMP2->GetYaxis()->SetLabelSize(0.035);
   PrevFitTMP2->GetYaxis()->SetTitleSize(0.035);
   PrevFitTMP2->GetYaxis()->SetTitleFont(42);
   PrevFitTMP2->SetParameter(0,0); // flat line at 1
   PrevFitTMP2->SetParError(0,0);
   PrevFitTMP2->SetParLimits(0,0,0);
   PrevFitTMP2->Draw("same");

  bottom->Modified();
  c1->cd();
  c1->Modified();
  c1->cd();

}

