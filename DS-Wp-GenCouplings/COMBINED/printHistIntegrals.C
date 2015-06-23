void printHistIntegrals() {

    TFile *infile = new TFile("/uscms_data/d2/dsperka/8TeV/MakeTBntuples/29Jan/CMSSW_5_3_6/src/UserCode/dsperka/wprimetb/RootFiles_For2DLimits_31May_finalbins_scaleGenTopPt/electron_BestJetJet2W_M_WprimeModRight_Histos-final_85_140_dr03_lep50.root","OLD");
  
  int i=0;

  TIter nextkey(gDirectory->GetListOfKeys());
  while (key = (TKey*)nextkey()) {
    obj = key->ReadObj(); 
    h = (TH1*)obj; 
    std::cout << h->GetName() << " " << h->Integral() <<  std::endl;

    ++i;
  }    

  std::cout<< i << " files" << std::endl;
}
