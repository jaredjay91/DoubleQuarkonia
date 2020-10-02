#include <iostream>
#include "TF2.h"
#include "TH2.h"
#include "TMath.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include <RooGlobalFunc.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"

double PDGmassJpsi = 3.0969;
double PDGmassUps = 9.46;

const float minjpsi = 2.0;
const float maxjpsi = 4.0;
const float minups = 6.0;
const float maxups = 14.0;

using namespace std;
using namespace RooFit;
void fitRooDataset() {

  gStyle->SetOptStat(0);

  TString inFileName = "dataset_DoubleMuon_ReReco_HLT13_4plmuons_opSign_Acc_pt35_SoftId_Vtx_Trig__other_opSign_sameAcc_SoftId_Vtx_20200921.root";
  TFile* inFile = TFile::Open(inFileName,"READ");

  TString RangeCut = Form("upsmass>%.2f && upsmass<%.2f && jpsimass>%.2f && jpsimass<%.2f",minups, maxups, minjpsi, maxjpsi );

  RooDataSet *dataset = (RooDataSet*)inFile->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("upsmass")), *(ws->var("jpsimass")) ), RangeCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  //delete dataset;
  cout << "####################################" << endl;
  ws->var("upsmass")->setRange(minups, maxups);
  ws->var("upsmass")->Print();
  ws->var("jpsimass")->setRange(minjpsi, maxjpsi);
  ws->var("jpsimass")->Print();

  int nUpsMassBin1d = 60;
  int nJpsiMassBin1d = 60;
  int nUpsMassBin2d = 13;
  int nJpsiMassBin2d = 13;
  TCanvas* c2 = new TCanvas("c2","c2",0,0,500,500);
  c2->cd();
  TH1* hh_data = ws->data("reducedDS")->createHistogram("upsmass,jpsimass",nUpsMassBin2d,nJpsiMassBin2d);
  hh_data->Draw("LEGO");

  
  //Build fit function for 2D fit

  //SIGNAL FUNCTIONS
  RooRealVar massUps("mUps","mean of the signal gaussian mass PDF",PDGmassUps, PDGmassUps-0.1, PDGmassUps+0.1 );
  RooRealVar sigmaUps("sigmaUps","width of the signal gaussian mass PDF in Upsilon axis",0.15, 0.05, 0.3);
  //RooGaussian* gaussUps = new RooGaussian("gaussUps","gaussian PDF for Upsilon", *(ws->var("upsmass")), massUps, sigmaUps);
  RooRealVar alphaCbUps("alphaCbUps","tail shift", 2.5, 1.0, 5.0);
  RooRealVar nCbUps("nCbUps","power order", 2.5, 1.0, 5.0);
  RooCBShape* cbUps = new RooCBShape("cbUps", "cystal Ball", *(ws->var("upsmass")), massUps, sigmaUps, alphaCbUps, nCbUps);
  //alphaCbUps.setVal(2.0);
  //alphaCbUps.setConstant();
  nCbUps.setVal(2.0);
  nCbUps.setConstant();

  RooRealVar massJpsi("mJpsi","mean of the signal gaussian mass PDF",PDGmassJpsi, PDGmassJpsi-0.05, PDGmassJpsi+0.05 );
  RooRealVar sigmaJpsi("sigmaJpsi","width of the signal gaussian mass PDF in Upsilon axis",0.03, 0.005, 0.1);
  //RooGaussian* gaussJpsi = new RooGaussian("gaussJpsi","gaussian PDF for Jpsi", *(ws->var("jpsimass")), massJpsi, sigmaJpsi);
  RooRealVar alphaCbJpsi("alphaCbJpsi","tail shift", 2.5, 1.0, 5.0);
  RooRealVar nCbJpsi("nCbJpsi","power order", 2.5, 1.0, 5.0);
  RooCBShape* cbJpsi = new RooCBShape("cbJpsi", "cystal Ball", *(ws->var("jpsimass")), massJpsi, sigmaJpsi, alphaCbJpsi, nCbJpsi);
  //alphaCbJpsi.setVal(2.0);
  //alphaCbJpsi.setConstant();
  nCbJpsi.setVal(1.5);
  nCbJpsi.setConstant();


  //BACKGROUND FUNCTIONS
  RooRealVar gaus_mu("gaus_mu","gaus_mu", 5, 3, 4) ;
  RooRealVar gaus_sigma("gaus_sigma","gaus_sigma", 1, 0.3, 2);
  RooGaussian* gaussBkg = new RooGaussian("gaussBkg","gaussian PDF in Jpsi axis", *(ws->var("jpsimass")), gaus_mu, gaus_sigma);

  RooRealVar err_mu("#mu","err_mu", 5, 0, 15) ;
  RooRealVar err_sigma("#sigma","err_sigma", 1, 0, 15);
  RooRealVar m_lambda("#lambda","m_lambda",  15, 0, 25);
  RooGenericPdf *erfBkg = new RooGenericPdf("erfBkg","Erf*Exp Background in Upsilon axis","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("upsmass")), m_lambda, err_mu, err_sigma) );


  //DOUBLE-QUARKONIA SIGNAL
  RooRealVar *nSig= new RooRealVar("nSig","double quarkonia signals",10,0,100);
  RooAbsPdf* sigPDF = new RooProdPdf("sigPDF", "double dimuon signal", *cbUps, *cbJpsi);


  //J/PSI + MUMU
  RooRealVar *nJpsiComb= new RooRealVar("nJpsiComb","Jpsi-combinatorial",100,0,1000);
  RooAbsPdf* jpsiComb = new RooProdPdf("jpsiComb", "jpsi+mumu signal", *cbJpsi, *erfBkg);

  //UPSILON + MUMU
  RooRealVar *nUpsComb= new RooRealVar("nUpsComb","Upsilon-combinatorial",100,0,1000);
  RooAbsPdf* upsComb = new RooProdPdf("upsComb", "ups+mumu signal", *cbUps, *gaussBkg);


  //PURE BACKGROUND
  RooRealVar *nBkg= new RooRealVar("nBkg","double dimuon background",1000,0,100000);
  RooAbsPdf* bkgPDF = new RooProdPdf("bkgPDF", "double dimuon background", *gaussBkg, *erfBkg);


  //BUILD THE MODEL
  RooAddPdf* model = new RooAddPdf( "model", "Sig + Bkg + JpsiComb + UpsComb", RooArgList(*sigPDF,*bkgPDF,*jpsiComb,*upsComb), RooArgList(*nSig,*nBkg,*nJpsiComb,*nUpsComb) );
  ws->import(*model);


  //Fit the model to the data
  RooFitResult* fitRes2 = (RooFitResult*)ws->pdf("model")->fitTo(*reducedDS,Save(),Hesse(kTRUE),Timer(kTRUE),Extended(kTRUE));
  ws->import(*fitRes2);

  //make plots
  TH1* hh_pdf = ws->pdf("model")->createHistogram("upsmass,jpsimass",40,40);
  hh_pdf->SetLineColor(kRed);
  cout << "Entries in hh_pdf = " << hh_pdf->GetEntries() << endl;
  //cout << "Entries in dataset = " << nSig->getVal()+nBkg->getVal() << endl;
  hh_pdf->Draw("same surf");

  TCanvas* c1 = new TCanvas("c1","c1",0,500,800,400);
  c1->Divide(2);

  c1->cd(1);
  RooPlot* myPlotJpsi = ws->var("jpsimass")->frame(nJpsiMassBin1d); // bins
  ws->data("reducedDS")->plotOn(myPlotJpsi);
  ws->pdf("model")->plotOn(myPlotJpsi);
  ws->pdf("model")->plotOn(myPlotJpsi,Name("bkgPDF"),Components(RooArgSet(*bkgPDF)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  ws->pdf("model")->plotOn(myPlotJpsi,Name("sigPDF"),Components(RooArgSet(*sigPDF)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlotJpsi,Name("jpsiComb"),Components(RooArgSet(*jpsiComb)),LineColor(kGreen+3),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlotJpsi,Name("upsComb"),Components(RooArgSet(*upsComb)),LineColor(6),LineWidth(2),LineStyle(2));
  myPlotJpsi->Draw();

  c1->cd(2);
  RooPlot* myPlotUps = ws->var("upsmass")->frame(nUpsMassBin1d); // bins
  ws->data("reducedDS")->plotOn(myPlotUps);
  ws->pdf("model")->plotOn(myPlotUps);
  ws->pdf("model")->plotOn(myPlotUps,Name("bkgPDF"),Components(RooArgSet(*bkgPDF)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  ws->pdf("model")->plotOn(myPlotUps,Name("sigPDF"),Components(RooArgSet(*sigPDF)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlotUps,Name("jpsiComb"),Components(RooArgSet(*jpsiComb)),LineColor(kGreen+3),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlotUps,Name("upsComb"),Components(RooArgSet(*upsComb)),LineColor(6),LineWidth(2),LineStyle(2));
  myPlotUps->Draw();

  TLegend* fitleg = new TLegend(0.6,0.6,0.89,0.89); fitleg->SetTextSize(12);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlotUps->findObject("reducedDS"),"Data","pe");
  //fitleg->AddEntry(myPlotUps->findObject("model"),"Total fit","l");
  fitleg->AddEntry(myPlotUps->findObject("sigPDF"),"J/#psi+#Upsilon Signal","l");
  fitleg->AddEntry(myPlotUps->findObject("jpsiComb"),"J/#psi+#mu^{+}#mu^{-}","l");
  fitleg->AddEntry(myPlotUps->findObject("upsComb"),"#Upsilon+#mu^{+}#mu^{-}","l");
  fitleg->AddEntry(myPlotUps->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");

  fitRes2->Print("v");

  double Ntot = reducedDS->sumEntries();
  double Njpsi = ws->var("nJpsiComb")->getVal();
  double Nups = ws->var("nUpsComb")->getVal();
  cout << "sumEntries = " << Ntot << endl;
  cout << "nJpsi = " << Njpsi << endl;
  cout << "nUps = " << Nups << endl;
  cout << "number of cross double-quarkonia events = " << Nups*Njpsi/Ntot << endl;

  c1->SaveAs("DoubleQuarkoniaFit_1DProjections.pdf");
  c1->SaveAs("DoubleQuarkoniaFit_1DProjections.png");
  c2->SaveAs("DoubleQuarkoniaFit_2D.pdf");
  c2->SaveAs("DoubleQuarkoniaFit_2D.png");

  TString outFileName = "DoubleQuarkoniaFitResults.root";
  TFile* outf = new TFile(outFileName,"recreate");
  c1->Write();
  c2->Write();
  ws->Write();
  outf->Close();

}
