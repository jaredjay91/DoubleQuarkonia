#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <iostream>

void makeDoubleQuarkoniaPlot(TString histoFileName="dimuonMassDoubleQuarkonia.root") {

  gStyle->SetOptStat(0);
  
  TChain* theChain = new TChain("myTree");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_PromptAOD_v1_Oniatree_addvn.root");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_PromptAOD_v2_Oniatree_addvn.root");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_4muons_PromptAOD_v1v2_Oniatree_addvn.root");
  theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_4muons_DoubleMuon_ReReco_Oniatree_addvn.root");

  // The histograms
  // sse variable size bins, multiplying the bin limit of a previous bin by a scale factor.
  // This will make the bins span about the same width in a log-x scale.
  
  cout << "dimuonYellowPlot: creating bin limit array..." << endl;

  //My log bins:
  double massmin = 0.0;
  double massmax = 100.0;
  double jpsimassmin = 2.0;
  double jpsimassmax = 4.0;
  double upsmassmin = 9.0;
  double upsmassmax = 10.5;
  double Zmassmin = 60.0;
  double Zmassmax = 120.0;
  const int nBins = 200;
  double bins[nBins+1];
  cout << "massbins = {";
  for (int im=0; im<nBins+1; im++) {
    bins[im] = pow(10,log10(massmax)*im/((double)nBins));
    cout << bins[im] << ", ";
  }
  cout << "};" << endl;

  cout << "dimuonYellowPlot: booking mass histogram with ";
  cout << nBins << " mass bins, mass limits: ";
  cout << bins[0] << " - " << bins[nBins] << endl;

  int kTrigSel = 13;
  TH1F* massHistoAllDimuons = new TH1F( "massHistoAllQuarkonia",  "massHistoAllQuarkonia;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", nBins, bins);
  TH1F* massHistoUpsilonCandidates = new TH1F( "massHistoUpsilonCandidates",  "massHistoUpsilonCandidates;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, upsmassmin, upsmassmax);
  TH1F* massHistoOtherDimuons = new TH1F( "massHistoOtherDimuons",  "massHistoOtherDimuons;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", nBins, bins);
  TH1F* massHistoOtherJpsis = new TH1F( "massHistoOtherJpsis",  "massHistoOtherJpsis;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, jpsimassmin, jpsimassmax);
  TH1F* massHistoOtherUpsilons = new TH1F( "massHistoOtherUpsilons",  "massHistoOtherUpsilons;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, upsmassmin, upsmassmax);
  TH1F* massHistoOtherZs = new TH1F( "massHistoOtherZs",  "massHistoOtherZs;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, Zmassmin, Zmassmax);

  cout << "dimuonYellowPlot: Draw histograms from tree..." << endl;

  int year = 2018; // choose year: 2011, 2013, 2015, 2018
  int triggerOption = 2; // 1= pp triggers, 2 = PbPb triggers

  Int_t nevent = theChain->GetEntries();

  Int_t           nTrig;
  Int_t           trigPrescale[26];   //[nTrig]
  ULong64_t       HLTriggers;
  Int_t           Reco_mu_size;
  Int_t           Reco_QQ_size;
  ULong64_t       Reco_QQ_trig[66];
  TClonesArray    *Reco_mu_4mom;
  TClonesArray    *Reco_QQ_4mom;
  Int_t           Reco_QQ_mupl_idx[1000];
  Int_t           Reco_QQ_mumi_idx[1000];
  Int_t           Centrality;
  Int_t           Reco_QQ_type[66];
  Int_t           Reco_QQ_sign[66];
  Float_t         Reco_QQ_VtxProb[1000];   //[Reco_QQ_size]
  Int_t           Reco_mu_nPixWMea[1000];
  Int_t           Reco_mu_nTrkWMea[1000];
  Float_t         Reco_mu_dxy[1000];
  Float_t         Reco_mu_dz[1000];

  TBranch *b_Reco_mu_size;
  TBranch *b_Reco_QQ_size;
  TBranch *b_Reco_mu_4mom;
  TBranch *b_Reco_QQ_4mom;
  TBranch *b_Reco_QQ_mupl_idx;
  TBranch *b_Reco_QQ_mumi_idx;
  TBranch *b_Reco_QQ_sign;
  TBranch *b_Centrality;
  TBranch *b_nTrig;
  TBranch *b_trigPrescale;
  TBranch *b_HLTriggers;
  TBranch *b_Reco_QQ_trig;
  TBranch *b_Reco_QQ_type;
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Reco_mu_nPixWMea;
  TBranch        *b_Reco_mu_nTrkWMea;
  TBranch        *b_Reco_mu_dxy;
  TBranch        *b_Reco_mu_dz;

  Reco_mu_4mom = 0;
  Reco_QQ_4mom = 0;

  theChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
  theChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
  theChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  theChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  theChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  theChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  theChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  theChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  theChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  theChain->SetBranchAddress("Reco_QQ_mupl_idx", &Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  theChain->SetBranchAddress("Reco_QQ_mumi_idx", &Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
  theChain->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
  theChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  theChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  theChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  theChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  theChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  theChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);

  Int_t           Reco_mu_SelectionType[1000];
  TBranch        *b_Reco_mu_SelectionType;
  theChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  int DIMUIDPASS = 0;
  int HLTPASS = 0;
  int VERTEXPASS = 0;
  int UPSCANDIDATES = 0;
  int OTHERDIMUONS = 0;

  int stepSize = nevent/1000;
  for(int i=0; i<nevent; i++){
    theChain->GetEvent(i);
    if(i%stepSize==0){cout<<">>>>> EVENT "<<i<<" / "<<theChain->GetEntries()<<" ("<<(int)(100.*i/theChain->GetEntries())<<"%)"<<endl;}
   //1. Find a high-pt-dimuon-triggered event with at least 4 reconstructed muons and at least 2 opposite-sign dimuons.
    for(int j=0; j<Reco_QQ_size; j++){

      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(j);
      TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[j]);
      TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[j]);

      if (!(Reco_QQ_type[j]==1 &&//100% of dimuons pass this "type"
          Reco_QQ_sign[j]==0 &&
          Reco_QQ_size>1 &&
          Reco_mu_size>3 &&
          RecoQQmupl->Pt()>3.5 &&
          RecoQQmumi->Pt()>3.5)){
        continue;
      }

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j]]&((int)pow(2,3)));

      //SoftMuon Cuts
      bool muplSoft = ( passMuonTypePl && 
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[j]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[j]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[j]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[j]])<20.) 
          ) ; 

      bool mumiSoft = ( passMuonTypeMi && 
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[j]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[j]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[j]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[j]])<20.)  
          ) ; 
      if ( !(muplSoft && mumiSoft) ) 
        continue; 

      DIMUIDPASS++;

      if( !( ((HLTriggers&     ((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel))) && ((Reco_QQ_trig[j]&((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel))) )
            ){
        continue;
      }
      HLTPASS++;

   //2. In each event, look for a really good upsilon candidate (dimuon with good vertex probability, in the upsilon mass range 8-14, etc.).
      if ( Reco_QQ_VtxProb[j]  < 0.01 ){
        continue;
      }
      VERTEXPASS++;

      massHistoAllDimuons->Fill(RecoQQ4mom->M());

      if (!(RecoQQ4mom->M()<upsmassmax && RecoQQ4mom->M()>upsmassmin)){
        continue;
      }
      UPSCANDIDATES++;
      massHistoUpsilonCandidates->Fill(RecoQQ4mom->M());

   //3. In each of those upsilon-candidate events, look for other dimuons among the other muons not related to the upsilon candidate.

      for(int j2=j+1; j2<Reco_QQ_size; j2++){
        TLorentzVector *RecoQQ4mom2 = (TLorentzVector*) Reco_QQ_4mom->At(j2);
        TLorentzVector *RecoQQmupl2 = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[j2]);
        TLorentzVector *RecoQQmumi2 = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[j2]);
        if (Reco_QQ_mupl_idx[j2]==Reco_QQ_mupl_idx[j] || Reco_QQ_mumi_idx[j2]==Reco_QQ_mumi_idx[j])
          continue;

        if (!(Reco_QQ_type[j2]==1 &&//100% of dimuons pass this "type"
            Reco_QQ_sign[j2]==0 //&&
            //RecoQQmupl2->Pt()>3.5 &&
            //RecoQQmumi2->Pt()>3.5
            ) ){
          continue;
        }

/*        if ( Reco_QQ_VtxProb[j2]  < 0.01 ){
          continue;
        }

        bool passMuonTypePl2 = true;
        passMuonTypePl2 = passMuonTypePl2 && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j2]]&((int)pow(2,1)));
        passMuonTypePl2 = passMuonTypePl2 && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j2]]&((int)pow(2,3)));

        bool passMuonTypeMi2 = true;
        passMuonTypeMi2 = passMuonTypeMi2 && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j2]]&((int)pow(2,1)));
        passMuonTypeMi2 = passMuonTypeMi2 && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j2]]&((int)pow(2,3)));

        //SoftMuon Cuts
        bool muplSoft2 = ( passMuonTypePl2 && 
            (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[j2]] > 5) &&
            (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[j2]] > 0) &&
            (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[j2]])<0.3) &&
            (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[j2]])<20.) 
            ) ; 

        bool mumiSoft2 = ( passMuonTypeMi2 && 
            (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[j2]] > 5) &&
            (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[j2]] > 0) &&
            (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[j2]])<0.3) &&
            (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[j2]])<20.)  
            ) ; 
        if ( !(muplSoft2 && mumiSoft2) ) 
          continue;

        //if( !( ((Reco_QQ_trig[j2]&((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel))) )
        //      ){
        //  continue;
        //}
*/
        OTHERDIMUONS++;
        massHistoOtherDimuons->Fill(RecoQQ4mom2->M());
        if (RecoQQ4mom2->M()<jpsimassmax && RecoQQ4mom2->M()>jpsimassmin){
          massHistoOtherJpsis->Fill(RecoQQ4mom2->M());
        }
        else if (RecoQQ4mom2->M()<upsmassmax && RecoQQ4mom2->M()>upsmassmin){
          massHistoOtherUpsilons->Fill(RecoQQ4mom2->M());
        }
        else if (RecoQQ4mom2->M()<Zmassmax && RecoQQ4mom2->M()>Zmassmin){
          massHistoOtherZs->Fill(RecoQQ4mom2->M());
        }
      }

   //4. Make 2 histograms of the upsilon candidates and the other good dimuons.

    }
  }

  //massHistoAllDimuons->Scale(1,"width");
  //massHistoUpsilonCandidates->Scale(1,"width");
  //massHistoOtherDimuons->Scale(1,"width");
  //massHistoOtherJpsis->Scale(1,"width");
  //massHistoOtherUpsilons->Scale(1,"width");
  //massHistoOtherZs->Scale(1,"width");

  cout << "DIMUIDPASS = " << DIMUIDPASS << endl;
  cout << "HLTPASS = " << HLTPASS << endl;
  cout << "VERTEXPASS = " << VERTEXPASS << endl;
  cout << "UPSCANDIDATES = " << UPSCANDIDATES << endl;
  cout << "OTHERDIMUONS = " << OTHERDIMUONS << endl;

  //make a quick plot
  TCanvas* c0 = new TCanvas("c0","c0",0,0,400,400);
  c0->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  massHistoAllDimuons->SetFillColor(kYellow);
  massHistoAllDimuons->SetTitle("dimuons");
  massHistoAllDimuons->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoAllDimuons->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoAllDimuons->GetXaxis()->CenterTitle();
  massHistoAllDimuons->Draw("hist");

  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,400);
  c1->Divide(2);
  c1->cd(1);
  massHistoUpsilonCandidates->SetFillColor(kYellow);
  massHistoUpsilonCandidates->SetTitle("Upsilon candidates");
  massHistoUpsilonCandidates->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoUpsilonCandidates->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoUpsilonCandidates->GetXaxis()->CenterTitle();
  massHistoUpsilonCandidates->Draw("hist");
  c1->cd(2);
  gPad->SetLogy();
  gPad->SetLogx();
  massHistoOtherDimuons->SetFillColor(kGreen+3);
  massHistoOtherDimuons->SetTitle("Other dimuons in UC events");
  massHistoOtherDimuons->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoOtherDimuons->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoOtherDimuons->GetXaxis()->CenterTitle();
  massHistoOtherDimuons->Draw("hist");

  //make a quick plot
  TCanvas* c2 = new TCanvas("c2","c2",0,0,800,800);
  c2->Divide(2,2);
  c2->cd(1);
  massHistoUpsilonCandidates->Draw("hist");
  c2->cd(2);
  massHistoOtherJpsis->SetFillColor(kGreen+3);
  massHistoOtherJpsis->SetTitle("Other Jpsis in UC events");
  massHistoOtherJpsis->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoOtherJpsis->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoOtherJpsis->GetXaxis()->CenterTitle();
  massHistoOtherJpsis->Draw("hist");
  c2->cd(3);
  massHistoOtherUpsilons->SetFillColor(kGreen+3);
  massHistoOtherUpsilons->SetTitle("Other Upsilons in UC events");
  massHistoOtherUpsilons->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoOtherUpsilons->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoOtherUpsilons->GetXaxis()->CenterTitle();
  massHistoOtherUpsilons->Draw("hist");
  c2->cd(4);
  massHistoOtherZs->SetFillColor(kGreen+3);
  massHistoOtherZs->SetTitle("Other Zs in UC events");
  massHistoOtherZs->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoOtherZs->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoOtherZs->GetXaxis()->CenterTitle();
  massHistoOtherZs->Draw("hist");
  TString cutsTag = "DoubleMuon_ReReco_HLT13_4muons_opSign_pt35_SoftId_Vtx_massCut__other_opSign_notBinScaled";
  c0->SaveAs(Form("Plots/DoubleQuarkonia0_%s.png",cutsTag.Data()));
  c0->SaveAs(Form("Plots/DoubleQuarkonia0_%s.pdf",cutsTag.Data()));
  c1->SaveAs(Form("Plots/DoubleQuarkonia1_%s.png",cutsTag.Data()));
  c1->SaveAs(Form("Plots/DoubleQuarkonia1_%s.pdf",cutsTag.Data()));
  c2->SaveAs(Form("Plots/DoubleQuarkonia2_%s.png",cutsTag.Data()));
  c2->SaveAs(Form("Plots/DoubleQuarkonia2_%s.pdf",cutsTag.Data()));

// Output file creation and writing of histogram
  cout << "dimuonDoubleQuarkoniaMakeHisto: Open output file... " << histoFileName << endl;
  TFile* outFile = new TFile(histoFileName,"RECREATE");

  cout << "Writing histograms .. " << endl;
  cout << massHistoUpsilonCandidates->GetName() << ", " << massHistoUpsilonCandidates->GetEntries() << " dimuons " << endl;
  massHistoUpsilonCandidates->Write();
  massHistoOtherDimuons->Write();
  massHistoOtherJpsis->Write();
  massHistoOtherUpsilons->Write();
  massHistoOtherZs->Write();

  outFile->Close();
  return;
}
