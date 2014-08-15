#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iomanip>

#include <TFile.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TCint.h>


 const int nmax=10;

TFile *vertex(Int_t get=0)
{

  TString filename = "vertex.root";
  TFile *hfile = 0;
  hfile = (TFile*)gROOT->FindObject(filename); if (hfile) hfile->Close();
  hfile = new TFile(filename,"RECREATE","Vertex Study");

  gBenchmark->Start("vertex");

  // Create a new canvas.                                                                      
  // Margins:
  gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadBottomMargin(0.13);
  gStyle -> SetPadLeftMargin(0.16);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetPadRightMargin(0.02);

  TCanvas *c1 = new TCanvas("c1","Vertex Plots",200,10,600,600);
  //c1->SetFillColor(42);
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  c1->SetFrameLineWidth(3);
  c1->SetRightMargin(0.05);


//  assume vertices are distributed as a gaussian with mean 0 and rms sigma
//  assume single vertices are never split
//  assume vertices that are closer than delta_z_merge are merged
// assume true PV multiplicity is a poisson with mean mean_num_int
//  assume there is 1 hard scatter and that the pileup is added on that
// with poisson distribution
  //Double_t delta_z_merge=1.5; // cm
  Double_t delta_z_merge=0.3; // cm
  Double_t sigma_z=6.1;  // cm
  //Double_t mean_num_int=0.6;
  //Double_t mean_num_int=0.49;
  //Double_t mean_num_int=1.56106;
  Double_t mean_num_int=1.45145;
  int nevents=13177;
  Double_t vtx[nmax];
  
 TRandom3 r;
 gRandom->SetSeed();

 TH1F *nvtxtrue  = new TH1F("nvtxtrue","true number of vertices",10,0.,10.);
 TH1F *ztrue = new TH1F("ztrue","true z distribution",100,-20.,20.);
 TH2F *nmt = new TH2F("nmt","ntrue vs nmeasured",10,0.,10.,10,0.,10.);

 for (int i=0;i<nevents;i++){  // make some events
   // pick number of vertices in this event
   Int_t nvtx = 1+r.Poisson(mean_num_int); // 1 is the h.s.  r is the pileup
   if(nvtx==2) std::cout<<std::endl;
   if(nvtx==2) std::cout<<"new event with true number of vertices "<<nvtx<<std::endl;
   float anvtx=nvtx-0.5;
   nvtxtrue->Fill(anvtx);
   // because I just don't program in C enough
   if(nvtx>=nmax) nvtx=nmax;
   // pick vertex position for each
   for(Int_t j=0;j<nvtx;j++) {
     vtx[j]=r.Gaus(0.,sigma_z);
     if(nvtx==2) std::cout<<"  vertex "<<j<<" at "<<vtx[j]<<std::endl;
     ztrue->Fill(vtx[j]);
   }
     // crude merge algorithm
   Int_t nvtx_measured = nvtx;
   if(nvtxtrue) {
   // start with last vertex and see if can merge with any earlier ones
     for(Int_t k=nvtx-1;k>0;k--) {
       int merge=0;
       for(Int_t l=0;l<k;l++) {
	 if(merge==0) {
           if(fabs(vtx[l]-vtx[k])<delta_z_merge) {
	     if(nvtx==2) std::cout<<" merging vertices "<<k<<" and "<<l<<std::endl;
	     merge++;
             vtx[l]=(vtx[k]+vtx[l])/2.;
	     for(Int_t m=0;m<nvtx_measured-k-1;m++){
               vtx[k+m]=vtx[k+m+1];
             }
	     nvtx_measured--;  
           }
         }
       }
     }
   }
   if(nvtx==2) std::cout<<"  number of measured vertices is "<<nvtx_measured<<std::endl;
   for(Int_t k=0;k<nvtx_measured;k++) {
     if(nvtx==2) std::cout<<"     measured vertex "<<k<<" at "<<vtx[k]<<std::endl;
   }
   float anvtx_measured = nvtx_measured-0.5;
   nmt->Fill(anvtx,anvtx_measured);
 }
 std::cout<<" summary of measured versus true "<<std::endl;
 float atest=0.;
 for(Int_t i=0;i<nmax;i++) {
   for(Int_t k=0;k<nmax;k++) {
     std::cout<<std::setw(8)<<std::setprecision(2)<<(nmt->GetBinContent(i,k))/nevents;
     atest+=(nmt->GetBinContent(i,k))/nevents;
   }
   std::cout<<std::endl;
 }
 std::cout<<"atest is "<<atest<<std::endl;

 std::cout<<" summary of events with 1 reconstructed vertex "<<std::endl;
 float asum=0.;
 for(Int_t i=0;i<nmax;i++) {
   asum+=nmt->GetBinContent(i,1);
 }
 std::cout<<" asum over nevets is "<<asum/nevents<<std::endl;
 for(Int_t i=0;i<nmax;i++) {
   std::cout<<"fraction of reconstructed 1 vertex events from true "<<i<<
     " vertex event is "<<nmt->GetBinContent(i,1)/asum<<std::endl;
 }


 gROOT->SetStyle("Plain");

 nvtxtrue->Draw();
 c1->Modified();
 c1->Update();

 gBenchmark->Show("vertex");

 // Save all objects in this file                                                             
 nvtxtrue->SetFillColor(0);
 hfile->Write();
 nvtxtrue->SetFillColor(400);
 c1->Modified();
 c1->SaveAs("vertices.pdf");
 return hfile;



}
