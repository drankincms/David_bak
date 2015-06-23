{
  //#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
  //#include<stdio.h>

using namespace std;



 TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR"); 
 // For the canvas:


 tdrStyle->SetCanvasBorderMode(0);
 tdrStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
 tdrStyle->SetCanvasDefH(600); //Height of canvas
 tdrStyle->SetCanvasDefW(600); //Width of canvas
 tdrStyle->SetCanvasDefX(0);   //POsition on screen
 tdrStyle->SetCanvasDefY(0);


 // For the Pad:
 tdrStyle->SetPadBorderMode(0);
 // tdrStyle->SetPadBorderSize(Width_t size = 1);
 tdrStyle->SetPadColor(0); // kWhite
 tdrStyle->SetPadGridX(0); //false
 tdrStyle->SetPadGridY(0); //false
 tdrStyle->SetGridColor(0);
 tdrStyle->SetGridStyle(3);
 tdrStyle->SetGridWidth(1);

 // For the frame:
 tdrStyle->SetFrameBorderMode(0);
 tdrStyle->SetFrameBorderSize(1);
 tdrStyle->SetFrameFillColor(0);
 tdrStyle->SetFrameFillStyle(0);
 tdrStyle->SetFrameLineColor(1);
 tdrStyle->SetFrameLineStyle(1);
 tdrStyle->SetFrameLineWidth(1);

 // For the histo:
 // tdrStyle->SetHistFillColor(1);
 // tdrStyle->SetHistFillStyle(0);
 tdrStyle->SetHistLineColor(1);
 tdrStyle->SetHistLineStyle(0);
 tdrStyle->SetHistLineWidth(1);
 // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
 // tdrStyle->SetNumberContours(Int_t number = 20);

 tdrStyle->SetEndErrorSize(2);
 //tdrStyle->SetErrorMarker(20);   /// I COMMENTED THIS OUT
 //tdrStyle->SetErrorX(0.);

 //tdrStyle->SetMarkerStyle(20);


 //For the fit/function:
 tdrStyle->SetOptFit(1011);
 tdrStyle->SetFitFormat("5.4g");
 tdrStyle->SetFuncColor(2);
 tdrStyle->SetFuncStyle(1);
 tdrStyle->SetFuncWidth(1);

 //For the date:
 tdrStyle->SetOptDate(0);
 // tdrStyle->SetDateX(Float_t x = 0.01);
 // tdrStyle->SetDateY(Float_t y = 0.01);

 // For the statistics box:
 tdrStyle->SetOptFile(0);
 tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
 tdrStyle->SetStatColor(0); // kWhite
 tdrStyle->SetStatFont(42);
 //tdrStyle->SetStatFontSize(0.025);
 tdrStyle->SetStatFontSize(0.04);
 tdrStyle->SetStatTextColor(1);
 tdrStyle->SetStatFormat("6.4g");
 tdrStyle->SetStatBorderSize(1);
 tdrStyle->SetStatH(0.1);
 tdrStyle->SetStatW(0.15);
 // tdrStyle->SetStatStyle(Style_t style = 1001);
 // tdrStyle->SetStatX(Float_t x = 0);
 // tdrStyle->SetStatY(Float_t y = 0);

 // For the Global title:

 tdrStyle->SetOptTitle(0);
 tdrStyle->SetTitleFont(42);
 tdrStyle->SetTitleColor(1);
 tdrStyle->SetTitleTextColor(1);
 tdrStyle->SetTitleFillColor(10);
 tdrStyle->SetTitleFontSize(0.05);
 // tdrStyle->SetTitleH(0); // Set the height of the title box
 // tdrStyle->SetTitleW(0); // Set the width of the title box
 // tdrStyle->SetTitleX(0); // Set the position of the title box
 // tdrStyle->SetTitleY(0.985); // Set the position of the title box
 // tdrStyle->SetTitleStyle(Style_t style = 1001);
 // tdrStyle->SetTitleBorderSize(2);

 // For the axis titles:

 tdrStyle->SetTitleColor(1, "XYZ");
 tdrStyle->SetTitleFont(42, "XYZ");
 tdrStyle->SetTitleSize(0.06, "XYZ");
 tdrStyle->SetTitleOffset(0.8, "XY"); // Another way to set the Offset
 tdrStyle->SetTitleOffset(1.27, "Z"); // Another way to set the Offset

 // For the axis labels:

 tdrStyle->SetLabelColor(1, "XYZ");
 tdrStyle->SetLabelFont(42, "XYZ");
 tdrStyle->SetLabelOffset(0.007, "XYZ");
 tdrStyle->SetLabelSize(0.04, "XYZ");

 // For the axis:

 tdrStyle->SetAxisColor(1, "XYZ");
 tdrStyle->SetStripDecimals(1); // kTRUE
 tdrStyle->SetTickLength(0.03, "XYZ");
 tdrStyle->SetNdivisions(510, "XY");
 tdrStyle->SetNdivisions(-510, "Z");
 tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
 tdrStyle->SetPadTickY(1);

 // Change for log plots:
 tdrStyle->SetOptLogx(0);
 tdrStyle->SetOptLogy(0);
 tdrStyle->SetOptLogz(0);

 // Postscript options:
 tdrStyle->SetPaperSize(20.,20.);





    
 /*
   tdrStyle->SetCanvasBorderMode(0);
   tdrStyle->SetCanvasColor(kWhite);
   tdrStyle->SetCanvasDefH(600); //Height of canvas
   tdrStyle->SetCanvasDefW(600); //Width of canvas
   tdrStyle->SetCanvasDefX(0);   //POsition on screen
   tdrStyle->SetCanvasDefY(0);

   // For the Pad:
   tdrStyle->SetPadBorderMode(0);
   // tdrStyle->SetPadBorderSize(Width_t size = 1);
   tdrStyle->SetPadColor(kWhite);
   tdrStyle->SetPadGridX(false);
   tdrStyle->SetPadGridY(false);
   tdrStyle->SetGridColor(0);
   tdrStyle->SetGridStyle(3);
   tdrStyle->SetGridWidth(1);

   // For the frame:
   tdrStyle->SetFrameBorderMode(0);
   tdrStyle->SetFrameBorderSize(1);
   tdrStyle->SetFrameFillColor(0);
   tdrStyle->SetFrameFillStyle(0);
   tdrStyle->SetFrameLineColor(1);
   tdrStyle->SetFrameLineStyle(1);
   tdrStyle->SetFrameLineWidth(1);

   // For the histo:
   // tdrStyle->SetHistFillColor(1);
   // tdrStyle->SetHistFillStyle(0);
   tdrStyle->SetHistLineColor(1);
   tdrStyle->SetHistLineStyle(0);
   tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("emr"); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // For the Global title:
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(22, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(22, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  //tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  //tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  //tdrStyle->SetPaintTextFormat(const char* format = "g");
  tdrStyle->SetPalette(1);
  //tdrStyle->SetTimeOffset(Double_t toffset);
  //tdrStyle->SetHistMinimumZero(kTRUE);
  */

    tdrStyle->SetPalette(1);
 
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  tdrStyle->SetNumberContours(NCont);

  //TLatex *lab = new TLatex(0.70,0.85, "CMS 2008");
  //lab->SetNDC();
  //lab->SetTextFont(42);
  //lab->SetTextSize(0.05);
  //lab->Draw("same");

  gROOT -> ForceStyle();

  tdrStyle->cd();





  Double_t m[14]={800,950,1100,1250,1400,1550,1700,1850,2000,2150,2300,2450,2500,2650};

  //  Double_t m[14]={800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400};

 TFile *file=new TFile("plots_expected.root");
 gStyle->SetTextFont(62);

 TH2F* obs1 = new TH2F(* ((TH2F*) file->Get("obs1")));

 // TH2 *obs1 = new TH2F("obs1","ar vs al", 101,-0.005,1.005, 101,-0.005, 1.005);

  
 TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
color=(TColor*)(gROOT->GetListOfColors()->At(100));
color->SetRGB(1,0.07,0.07);
color=(TColor*)(gROOT->GetListOfColors()->At(99));
color->SetRGB(1,0.14,0.14);
color=(TColor*)(gROOT->GetListOfColors()->At(98));
color->SetRGB(1,0.21,0.21);
color=(TColor*)(gROOT->GetListOfColors()->At(97));
color->SetRGB(1,0.28,0.28);
color=(TColor*)(gROOT->GetListOfColors()->At(96));
color->SetRGB(1,0.35,0.35);
color=(TColor*)(gROOT->GetListOfColors()->At(95));
color->SetRGB(1,0.42,0.42);
color=(TColor*)(gROOT->GetListOfColors()->At(94));
color->SetRGB(1,0.49,0.49);
color=(TColor*)(gROOT->GetListOfColors()->At(93));
color->SetRGB(1,0.56,0.56);
color=(TColor*)(gROOT->GetListOfColors()->At(92));
color->SetRGB(1,0.63,0.63);
color=(TColor*)(gROOT->GetListOfColors()->At(91));
color->SetRGB(1,0.70,0.70);
color=(TColor*)(gROOT->GetListOfColors()->At(90));
color->SetRGB(1,0.79,0.79);
Int_t color_indices[11]={90,91,92,93,94,95,96,97,98,99,100};
gStyle->SetPalette(11,color_indices);

 c1->SetFillColor(kWhite);                                                                                           
 obs1->Draw("cont3zcol");
 obs1->SetLineWidth(3);
 c1->SetTopMargin(0.1);
 c1->SetBottomMargin(0.12);
 c1->SetRightMargin(0.2);
 c1->SetLeftMargin(0.12);

 obs1->SetContour(11,m);
 // obs1->Draw("cont2same");  

 obs1->GetXaxis()->SetTitle("a^{L}");

 obs1->GetYaxis()->SetTitle("a^{R}");

 obs1->GetZaxis()->SetTitle("M(W') [GeV]");
 obs1->GetZaxis()->SetRangeUser(800,2300);

TLatex *t1 = new TLatex();
float x1=0.2;
float y1=.85;

t1->SetNDC();
t1->SetTextFont(42);
t1->SetTextColor(1);
t1->SetTextAlign(11);
t1->SetTextSize(0.06);
t1->SetTextFont(42);

//TLatex *t1temp=0;

 TLatex *latex2=new TLatex();
 latex2->SetNDC();
 latex2->SetTextSize(0.04);
 latex2->SetTextFont(62);
 latex2->SetTextAlign(31); // align right
 //latex2->DrawLatex(0.802, 0.914, "CMS Preliminary, 19.6 fb^{-1} at #sqrt{s} = 8 TeV");
 latex2->DrawLatex(0.802, 0.914, "CMS, L=19.5 fb^{-1} at #sqrt{s} = 8 TeV");
 latex2->SetTextAlign(11); // align left
 //latex2->DrawLatex(0.12, 0.92, "CMS Simulation");


 TLatex *latex1=new TLatex();
 latex1->SetNDC();
 latex1->SetTextSize(0.04);
 latex1->SetTextAlign(31); // align right
 // latex1->DrawLatex(0.85, 0.86, " M(W') [GeV]");


 TLatex *latex4=new TLatex();
 latex4->SetNDC();
 latex4->SetTextSize(0.045);
 latex4->SetTextFont(62);
 latex4->SetTextAlign(22); // center
 //latex2->DrawLatex(0.53, 0.25, "e/#mu+jets 95% C.L. Expected");
 latex4->DrawLatex(0.34, 0.19, "e/#mu+jets N_{b tags} = 1 or 2");
 latex4->DrawLatex(0.34, 0.25, "95% CL expected");



c1->Update();
c1->Print("contour_arVSal-combined_expected150_final.pdf");
//c1->Print("contour_arVSal-combined_expected150.C");
c1->Print("contour_arVSal-combined_expected150_final.png");

    
}//file
  
