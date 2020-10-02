 
void skimEvents2015()
{

  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262620_263757_noCUT.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0B_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root";
  TString fnameData3 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0C_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root";
  TString fnameData4 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0D_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root";

  // Get old file, old tree and set top branch address
  TChain* theChain = new TChain("hionia/myTree");
  theChain->Add(fnameData1.Data());
  theChain->Add(fnameData2.Data());
  theChain->Add(fnameData3.Data());
  theChain->Add(fnameData4.Data());

  //TFile* file = TFile::Open("PromptAOD_v1_Oniatree_addvn_part2.root");
  //TTree* theChain = (TTree*)file->Get("myTree");

  int kTrigSel = 0;
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

  const int nentries = theChain->GetEntries();
  //const int nentries = 100000;
  cout << "nentries = " << nentries << endl;
  const int stepSize = nentries/1000;

  //Event *event = nullptr;
  //theChain->SetBranchAddress("event", &event);
 
  // Create a new file + a clone of old tree in new file
  TFile newfile("skimmed_HLT0_4muons_DoubleMuon_PromptReco_Oniatree2015.root", "recreate");
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
