//
// PileUp Standalone
//
// Zaixing Mao  
// 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "pileupFly.h"

//pileupSkim works as a standalone calculator with no needs for the CMS frame work
//MC and Data distributions are stored in "pileupFly.h" file under the USER OPTIONS section
//pileupSkim takes in the above distributions and the number of true interactions to calcuate the pileup reweight 
//The new PU reweight is stored in a new root file under branch "weight_PileUp_new" by default

//One can also do the reweighting on the fly in your own event loop by: Include pileupFly.h and copy the two lines highlighted by @@@@@

int pileupSkim() {

    
  //**************************************************
  //*******            USER OPTIONS          *********
  //**************************************************
  bool const newFile = true;                          //The Default Saving Choice is to Create a new file with _new following the old name
                                                      //Set false to overwrite input file. File ofileName.root will still be created but empty.
  
  TString ifileName = "ljmet_tree";                   //Input file name
  
  TString ofileName = ifileName + "_new";             //Output file name
  
  TString nTrueName = "nTrueInteractions_PileUpCalc"; //Name of the branch storing nTrueInteractions value in the input tree.
    
  TString newBranchName = "weight_PileUp_new";        //New branch for storing pileup weightings
  
  bool const reNameBranch = true;                     //If found branch with name newBranchName in input file:
                                                      //true (Default): create branch with name newBranchName + "_new" instead.
                                                      //false: abort process.
  //**************************************************
    
  //___Define Variables___
  double nTrue = -1;
  double MyWeight;
  TBranch* weightPU;
  TTree* otree;
  TString readOption = "READ";
  
    
  //@@@@@___Set MC and Data Distribution___@@@@@
  PU::SetDist();
    
    
  //___Check If Overwritting the Input File___
  if(!newFile) {
     readOption = "update";
     std::cout<<"Warning :: Overwritting input file  \""<<ifileName<<".root\"!!!"<<std::endl;
  }
   
    
  //___Input Tree___
  TFile* ifile = new TFile(ifileName + ".root",readOption);
  TTree* itree = (TTree*)ifile->Get("ljmet");
  
    
  //___Check If Branch "newBranchName" Already Exist in Input File___
  while (itree->FindBranch(newBranchName)!=NULL){
     std::cout<<"Warning :: Branch :\""<<newBranchName<<"\" already exist in file \""<<ifileName<<".root\"!!!"<<std::endl;
     
     //_____Check If Abort_____
     if(!reNameBranch) {std::cout<<"Abort"<<std::endl; return 0;}
     
     newBranchName+="_new";
     std::cout<<"Try     :: Branch :\""<<newBranchName<<"\" "<<std::endl;
  }
    
    
  //___Output Tree___
  std::cout<<"Creating:: File: \""<<ofileName<<".root"<<"\" ..." <<std::endl;
  TFile* ofile = new TFile(ofileName+ ".root","RECREATE");
  if(newFile) otree = (TTree*)itree->CloneTree();
    
    
  //___Create Weight PU Branch___
  std::cout<<"Booking :: Branch: \""<<newBranchName<<"\" ..." <<std::endl;
  if(newFile) weightPU = otree->Branch(newBranchName,&MyWeight,newBranchName + "/D");
  if(!newFile) weightPU = itree->Branch(newBranchName,&MyWeight,newBranchName + "/D");
  
    
  //___Get Necessary Branches___
  itree->SetBranchAddress(nTrueName,&nTrue);
  

  //___Event Loop___
  int nEvent= itree->GetEntries();
  for(int ev = 0; ev<nEvent; ev++){
    itree->GetEvent(ev);
    
    //@@@@@___Calculate LumiReWeights___@@@@@
    MyWeight = PU::LumiReWeights(nTrue);    
    
    weightPU->Fill();
  }
  
  
  //___Write Files___
  if(newFile) ofile->Write();
  if(!newFile) {
    itree->Write("", TObject::kOverwrite);
    ifile->Write("",TObject::kOverwrite);
  }
  delete ifile;
  delete ofile;
  return 0;
}

