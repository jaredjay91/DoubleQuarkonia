 
void skimEvents()
{

  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameDataReReco = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part*.root";

  // Get old file, old tree and set top branch address
  TChain* theChain = new TChain("myTree");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/PromptAOD_v1_Oniatree_addvn_part1.root");
  //theChain->Add("/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbRD_Upsilonv2/PromptAOD_v1_Oniatree_addvn_part2.root");
  //theChain->Add(fnameDataReReco.Data()); 
  //theChain->Add(fnameData1.Data());
  //theChain->Add(fnameData2.Data());
  theChain->Add(fnameDataReReco.Data());

  //TFile* file = TFile::Open("PromptAOD_v1_Oniatree_addvn_part2.root");
  //TTree* theChain = (TTree*)file->Get("myTree");

  int kTrigSel = 13;
  const int maxBranchSize = 1000;

  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       HLTriggers;
  Int_t           Reco_mu_size;

  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_mu_size;

  theChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  theChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  theChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);

  //const int nentries = theChain->GetEntries();
  const int nentries = 1000000;
  cout << "nentries = " << nentries << endl;
  const int stepSize = nentries/1000;

  //Event *event = nullptr;
  //theChain->SetBranchAddress("event", &event);
 
  // Create a new file + a clone of old tree in new file
  TFile newfile("skimmed_10millionevents_4muons_HLT13_DoubleMuon_ReReco_Oniatree_addvn.root", "recreate");
  //TFile newfile("/eos/cms/store/group/phys_heavyions/dileptons/jjay/DoubleQuarkonia/skimmedHLT13_4muons_PromptAOD_v1v2_Oniatree_addvn.root", "recreate");
  auto newtree = theChain->CloneTree(0);
 
  cout << "Now entering loop." << endl;
  for(int iev=0; iev<nentries; ++iev)
  {
    if(iev%stepSize==0) cout << ">>>>> EVENT " << iev << " / " << theChain->GetEntries() <<  " ("<<(int)(100.*iev/nentries) << "%)" << endl;
    //cout << "Cloning entry " << iev << endl;
    theChain->GetEntry(iev);

    if(Reco_mu_size<4) continue;
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    //if(!( (Reco_QQ_trig[iev]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;//this is a dimuon trigger, not an event trigger.

    newtree->Fill();
    //event->Clear();
  }

  const int newentries = newtree->GetEntries();
  cout << "new tree has " << newentries << " events" << endl;
  cout << "That's a " << 100.0*(nentries-newentries)/(double)nentries  << "% reduction!" << endl;
  //newtree->Print();
  newfile.Write();
}
