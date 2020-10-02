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


bool AcceptanceCutLoose(TLorentzVector* Muon) {
  return (
  (fabs(Muon->Eta())<2.4) && (
  (fabs(Muon->Eta())<0.3 && Muon->Pt()>3.4) ||
  (0.3<=fabs(Muon->Eta()) && fabs(Muon->Eta())<1.1 && Muon->Pt()>3.3) ||
  (1.1<=fabs(Muon->Eta()) && fabs(Muon->Eta())<1.4 && Muon->Pt()>7.7-4.0*fabs(Muon->Eta())) ||
  (1.4<=fabs(Muon->Eta()) && fabs(Muon->Eta())<1.55 && Muon->Pt()>2.1) ||
  (1.55<=fabs(Muon->Eta()) && fabs(Muon->Eta())<2.2 && Muon->Pt()>4.25-1.39*fabs(Muon->Eta())) ||
  (2.2<=fabs(Muon->Eta()) && Muon->Pt()>1.2)
  ));
};

bool AcceptanceCutTight(TLorentzVector* Muon) {
  return (
  (fabs(Muon->Eta())<2.4) && (
  (fabs(Muon->Eta())<1.2 && Muon->Pt()>3.5) ||
  (1.2<=fabs(Muon->Eta()) && fabs(Muon->Eta())<2.1 && Muon->Pt()>5.47-1.89*fabs(Muon->Eta())) ||
  (2.1<=fabs(Muon->Eta()) && Muon->Pt()>1.5)
  ));
};


void make2DmassPlot(int dateStr=20200921) {

  gStyle->SetOptStat(0);
  
  TChain* theChain = new TChain("myTree");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_PromptAOD_v1_Oniatree_addvn.root");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_PromptAOD_v2_Oniatree_addvn.root");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_4muons_PromptAOD_v1v2_Oniatree_addvn.root");
  theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/skimmedHLT13_4muons_DoubleMuon_ReReco_Oniatree_addvn.root");

  // The histogram
  double jpsimassmin = 2.0;
  double jpsimassmax = 4.0;
  double upsmassmin = 6.0;
  double upsmassmax = 14;
  const int njpsibins = 12;
  const int nupsbins = 13;

  double tightjpsimassmin = 3.0;
  double tightjpsimassmax = 3.2;
  double tightupsmassmin = 9.0;
  double tightupsmassmax = 10.0;

  TH2F* massHisto2D = new TH2F( "massHisto2D",  "massHisto2D", nupsbins, upsmassmin, upsmassmax, njpsibins, jpsimassmin, jpsimassmax);
  TH1F* massHistoUpsilons = new TH1F( "massHistoUpsilons",  "massHistoUpsilons;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, upsmassmin, upsmassmax);
  TH1F* massHistoJpsis = new TH1F( "massHistoJpsis",  "massHistoJpsis;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, jpsimassmin, jpsimassmax);
  TH1F* massHistoUpsilonsWithTightJpsi = new TH1F( "massHistoUpsilonsWithTightJpsi",  "massHistoUpsilonsWithTightJpsi;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, upsmassmin, upsmassmax);
  TH1F* massHistoJpsisWithTightUpsilon = new TH1F( "massHistoJpsisWithTightUpsilon",  "massHistoJpsisWithTightUpsilon;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, jpsimassmin, jpsimassmax);

  int kTrigSel = 13;

  Int_t           nTrig;
  UInt_t           eventNb;
  UInt_t           runNb;
  UInt_t           LS;
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
  Int_t           Ntracks;
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
  TBranch *b_Ntracks;
  TBranch *b_nTrig;
  TBranch *b_eventNb;
  TBranch *b_runNb;
  TBranch *b_LS;
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
  theChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  theChain->SetBranchAddress("runNb", &runNb, &b_runNb);
  theChain->SetBranchAddress("LS", &LS, &b_LS);
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
  theChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
  theChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  theChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  theChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  theChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  theChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);

  Int_t           Reco_mu_SelectionType[1000];
  TBranch        *b_Reco_mu_SelectionType;
  theChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Bool_t          Reco_mu_highPurity[1000];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  theChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Bool_t          Reco_mu_TrkMuArb[1000];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_TrkMuArb;   //!
  theChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);

  Bool_t          Reco_mu_TMOneStaTight[1000];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  theChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);

  TString cutsTag = "DoubleMuon_ReReco_HLT13_4plmuons_opSign_Acc_pt35_SoftId_Vtx_Trig__other_opSign_sameAcc_SoftId_Vtx";
  TFile* newfile;
  newfile = new TFile(Form("dataset_%s_%i.root",cutsTag.Data(),dateStr),"recreate");

  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* upsmassVar  = new RooRealVar("upsmass","upsilon mass variable",upsmassmin,upsmassmax,"GeV/c^{2}");
  RooRealVar* jpsimassVar  = new RooRealVar("jpsimass","jpsi mass variable",jpsimassmin,jpsimassmax,"GeV/c^{2}");
  RooArgSet* argSet    = new RooArgSet(*upsmassVar, *jpsimassVar);
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  int DIMUIDPASS = 0;
  int HLTPASS = 0;
  int VERTEXPASS = 0;
  int UPSCANDIDATES = 0;
  int JPSICANDIDATES = 0;
  int DIDIMUONS = 0;

  newfile->cd();

  Int_t nevent = theChain->GetEntries();
  //Int_t nevent = 100000;
  int stepSize = nevent/100;
  for(int i=0; i<nevent; i++){
    theChain->GetEvent(i);
    //if(i%stepSize==0){cout<<">>>>> EVENT "<<i<<" / "<<theChain->GetEntries()<<" ("<<(int)(100.*i/theChain->GetEntries())<<"%)"<<endl;}
   //1. Find a high-pt-dimuon-triggered event with at least 4 reconstructed muons and at least 2 opposite-sign dimuons.

    if(!(Reco_QQ_size>1 && Reco_mu_size>3) ){
      continue;
    }

    if( !((HLTriggers&     ((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel)) )
          ){
      continue;
    }

    for(int j=0; j<Reco_QQ_size; j++){
      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(j);
      TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[j]);
      TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[j]);

      if (!(abs(RecoQQ4mom->Rapidity())<2.4 &&
          RecoQQ4mom->Pt()<50
        )){
        continue;
      }
      if (!(//Reco_QQ_type[j]==1 &&//100% of dimuons pass this "type"
          Reco_QQ_sign[j]==0 &&
          //abs(Reco_QQ_sign[j])>0 &&
          //Reco_mu_size==4 &&
          RecoQQmupl->Pt()>3.5 &&
          RecoQQmumi->Pt()>3.5 /*&&
          abs(RecoQQmupl->Eta())<2.4 &&
          abs(RecoQQmumi->Eta())<2.4*/)){
        continue;
      }
      if (!( AcceptanceCutLoose(RecoQQmumi) && AcceptanceCutLoose(RecoQQmupl) )) {
        continue;
      }

      bool passMuonTypePl = true;
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j]]&((int)pow(2,1)));
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j]]&((int)pow(2,1)));
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j]]&((int)pow(2,3)));

      //SoftMuon Cuts
      bool muplSoft = ( passMuonTypePl && 
          //Reco_mu_highPurity[Reco_QQ_mupl_idx[j]] &&
          //Reco_mu_TrkMuArb[Reco_QQ_mupl_idx[j]] &&
          //Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[j]] &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[j]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[j]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[j]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[j]])<20.) 
          ) ; 

      bool mumiSoft = ( passMuonTypeMi && 
          //Reco_mu_highPurity[Reco_QQ_mumi_idx[j]] &&
          //Reco_mu_TrkMuArb[Reco_QQ_mumi_idx[j]] &&
          //Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[j]] &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[j]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[j]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[j]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[j]])<20.)  
          ) ; 
      if ( !(muplSoft && mumiSoft) ) 
        continue; 

      DIMUIDPASS++;

      if( !((Reco_QQ_trig[j]&((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel)))  ){
        continue;
      }
      HLTPASS++;

   //2. In each event, look for a really good upsilon candidate (dimuon with good vertex probability, in the upsilon mass range 8-14, etc.).
      if ( Reco_QQ_VtxProb[j] < 0.01 ){
        continue;
      }
      VERTEXPASS++;

      if (!(RecoQQ4mom->M()<upsmassmax && RecoQQ4mom->M()>upsmassmin)){
        continue;
      }
      UPSCANDIDATES++;
      //massHistoUpsilons->Fill(RecoQQ4mom->M());

   //3. In each of those upsilon-candidate events, look for other dimuons among the other muons not related to the upsilon candidate.

      for(int j2=0; j2<Reco_QQ_size; j2++){
        if (j2==j) continue;

        TLorentzVector *RecoQQ4mom2 = (TLorentzVector*) Reco_QQ_4mom->At(j2);
        TLorentzVector *RecoQQmupl2 = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[j2]);
        TLorentzVector *RecoQQmumi2 = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[j2]);

if( ((Reco_QQ_trig[j2]&((ULong64_t)pow(2,kTrigSel)))==((ULong64_t)pow(2,kTrigSel)) )
              ){
          //continue;
        }

        if (Reco_QQ_mupl_idx[j2]==Reco_QQ_mupl_idx[j] || Reco_QQ_mumi_idx[j2]==Reco_QQ_mumi_idx[j])
          continue;
        if (!(abs(RecoQQ4mom2->Rapidity())<2.4 &&
            RecoQQ4mom2->Pt()<50
          )){
          continue;
        }
        if (!(//Reco_QQ_type[j2]==1 &&//100% of dimuons pass this "type"
            Reco_QQ_sign[j2]==0 /*&&
            RecoQQmupl2->Pt()>3.5 &&
            RecoQQmumi2->Pt()>3.5 &&
            abs(RecoQQmupl2->Eta())<2.4 &&
            abs(RecoQQmumi2->Eta())<2.4*/)){
          continue;
        }
        if (!( AcceptanceCutLoose(RecoQQmumi2) && AcceptanceCutLoose(RecoQQmupl2) )) {
          continue;
        }

        if ( Reco_QQ_VtxProb[j2] < 0.01 ){
          continue;
        }

        bool passMuonTypePl2 = true;
        //passMuonTypePl2 = passMuonTypePl2 && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j2]]&((int)pow(2,1)));
        //passMuonTypePl2 = passMuonTypePl2 && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[j2]]&((int)pow(2,3)));

        bool passMuonTypeMi2 = true;
        //passMuonTypeMi2 = passMuonTypeMi2 && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j2]]&((int)pow(2,1)));
        //passMuonTypeMi2 = passMuonTypeMi2 && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[j2]]&((int)pow(2,3)));

        //SoftMuon Cuts
        bool muplSoft2 = ( passMuonTypePl2 && 
            //Reco_mu_highPurity[Reco_QQ_mupl_idx[j2]] &&
            //Reco_mu_TrkMuArb[Reco_QQ_mupl_idx[j2]] &&
            //Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[j2]] &&
            (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[j2]] > 5) &&
            (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[j2]] > 0) &&
            (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[j2]])<0.3) &&
            (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[j2]])<20.) 
            ) ; 

        bool mumiSoft2 = ( passMuonTypeMi2 && 
            //Reco_mu_highPurity[Reco_QQ_mumi_idx[j2]] &&
            //Reco_mu_TrkMuArb[Reco_QQ_mumi_idx[j2]] &&
            //Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[j2]] &&
            (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[j2]] > 5) &&
            (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[j2]] > 0) &&
            (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[j2]])<0.3) &&
            (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[j2]])<20.)  
            ) ; 
        if ( !(muplSoft2 && mumiSoft2) ) 
          continue;

        if (RecoQQ4mom2->M()<jpsimassmax && RecoQQ4mom2->M()>jpsimassmin){
          JPSICANDIDATES++;
          massHistoJpsis->Fill(RecoQQ4mom2->M());
          massHistoUpsilons->Fill(RecoQQ4mom->M());
          DIDIMUONS++;
          massHisto2D->Fill(RecoQQ4mom->M(),RecoQQ4mom2->M());
          //Print information on the dimuons going into the plot:
          cout << runNb << ":" << LS << ":" << eventNb << ", ";
          cout << "Cent " << Centrality << ", ";
          cout << "NT " << Ntracks << ", < ";
          cout << Form("m%.2f",RecoQQ4mom->M()) << ",";
          cout << "(" << Reco_QQ_mupl_idx[j] << "," << Reco_QQ_mumi_idx[j] << "); ";
          cout << Form("m%.2f",RecoQQ4mom2->M()) << ",";
          cout << "(" << Reco_QQ_mupl_idx[j2] << "," << Reco_QQ_mumi_idx[j2] << "); ";
          cout << ">" << endl;
          
          //Fill tight-mass histograms for 1D fitting:
          if (RecoQQ4mom2->M()<tightjpsimassmax && RecoQQ4mom2->M()>tightjpsimassmin){
            massHistoUpsilonsWithTightJpsi->Fill(RecoQQ4mom->M());
          }
          if (RecoQQ4mom->M()<tightupsmassmax && RecoQQ4mom->M()>tightupsmassmin){
            massHistoJpsisWithTightUpsilon->Fill(RecoQQ4mom2->M());
          }

          upsmassVar->setVal( (double)RecoQQ4mom->M() ) ;
          jpsimassVar->setVal( (double)RecoQQ4mom2->M() ) ;
          dataSet->add( *argSet);
        }
      }

    }
  }

  cout << "DIMUIDPASS = " << DIMUIDPASS << endl;
  cout << "HLTPASS = " << HLTPASS << endl;
  cout << "VERTEXPASS = " << VERTEXPASS << endl;
  cout << "UPSCANDIDATES = " << UPSCANDIDATES << endl;
  cout << "JPSICANDIDATES = " << JPSICANDIDATES << endl;
  cout << "DIDIMUONS = " << DIDIMUONS << endl;

  //make a quick plot
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,400);
  c1->Divide(2);
  c1->cd(1);
  massHistoUpsilons->SetFillColor(kYellow);
  massHistoUpsilons->SetTitle("Upsilon candidates");
  massHistoUpsilons->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoUpsilons->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoUpsilons->GetXaxis()->CenterTitle();
  massHistoUpsilons->Draw("hist");
  c1->cd(2);
  massHistoJpsis->SetFillColor(kGreen+3);
  massHistoJpsis->SetTitle("Jpsi candidates in UC events");
  massHistoJpsis->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoJpsis->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoJpsis->GetXaxis()->CenterTitle();
  massHistoJpsis->Draw("hist");

  //make a quick plot
  TCanvas* c2 = new TCanvas("c2","c2",500,0,500,500);
  massHisto2D->Draw("LEGO2");
  TCanvas* c3 = new TCanvas("c3","c3",0,500,500,500);
  massHisto2D->Draw("BOX");
  TCanvas* c4 = new TCanvas("c4","c4",500,500,500,500);
  massHisto2D->Draw("COLZ");

  c1->SaveAs(Form("Plots/2DmassPlotCandidates_%s.png",cutsTag.Data()));
  c1->SaveAs(Form("Plots/2DmassPlotCandidates_%s.pdf",cutsTag.Data()));
  c2->SaveAs(Form("Plots/2DmassPlot%s.png",cutsTag.Data()));
  c2->SaveAs(Form("Plots/2DmassPlot%s.pdf",cutsTag.Data()));
  c3->SaveAs(Form("Plots/2DmassBoxPlot%s.png",cutsTag.Data()));
  c3->SaveAs(Form("Plots/2DmassBoxPlot%s.pdf",cutsTag.Data()));
  c4->SaveAs(Form("Plots/2DmassColzPlot%s.png",cutsTag.Data()));
  c4->SaveAs(Form("Plots/2DmassColzPlot%s.pdf",cutsTag.Data()));

// Output file creation and writing of histogram
  newfile->cd();

  cout << "Writing histograms .. " << endl;
  massHistoUpsilons->Write();
  massHistoJpsis->Write();
  massHistoUpsilonsWithTightJpsi->Write();
  massHistoJpsisWithTightUpsilon->Write();
  massHisto2D->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();

  dataSet->Write(); 
  newfile->Close();

  return;
}
