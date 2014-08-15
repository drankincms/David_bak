#include<string.h>
using std::string;

void plotLimits(string limitDIR,  string limitFILE,  string titchan)
{

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
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

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.15);
  tdrStyle->SetPadRightMargin(0.05);

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
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
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

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  tdrStyle->SetNumberContours(NCont);

  gROOT -> ForceStyle();

  tdrStyle->cd();

  float xxx,yyy,charsize=0.07;

  gStyle->SetErrorX(5); 
  gStyle->SetEndErrorSize(5); 

  float X[22],EXP[22],EXPA[22],EXPB[22],EXPM1[22],EXPP1[22],EXPM2[22],EXPP2[22],OBS[22],CS[22],CS2[22];

  //  char name[10],expfile[80],obsfile[80],theoryfile[80],dirname[10];
  string expfile, expfileA, expfileB, obsfile,theoryfile, limitName, thTYPE;

  cout << " file name " << limitFILE << endl;

  //string explim = "_cls_limits_expected.txt";
  //string obslim = "_cls_limits_observed.txt";  
  string explim = "_bayesian_limits_expected.txt";
  string obslim = "_bayesian_limits_observed.txt";
  string thcurve = "_theory.txt";
  string sRH="_R_";
  string sLH="_L_";
  string sMix="_RL_";
  string mu="mu";
  string elec="elec";

  expfile=limitDIR+limitFILE+explim;
  obsfile=limitDIR+limitFILE+obslim;

  cout << expfile << endl;
  cout << obsfile << endl;

  //if (limitFILE.find(sRH)<30) {thTYPE="R";}
  //if (limitFILE.find(sLH)<30) {thTYPE="L";}
  //if (limitFILE.find(sMix)<30) {thTYPE="LR";}

  thTYPE="R";

  if (limitFILE.find(mu)>=40 && limitFILE.find(elec)>=40) const char *leptonType = "e/#mu+jets N_{b tags} #geq 1";
  if (limitFILE.find(mu)<40) const char *leptonType="#mu+jets N_{b tags} #geq 1";
  if (limitFILE.find(elec)<40) const char *leptonType="e+jets N_{b tags} #geq 1";
  theoryfile=thTYPE+thcurve;
  cout << " th file name " << theoryfile << endl;
  cout << " exp file name " << expfile << endl;
  cout << " obs file name " << obsfile << endl;

  int i=0;
  ifstream expFile;
  expFile.open(expfile.c_str());

  i=0;
  double xmass=0.;
  double maxmass = 3000.;
  while( !expFile.eof() ) {
    expFile >> X[i] >> EXP[i] >> EXPM2[i] >> EXPP2[i] >> EXPM1[i] >> EXPP1[i];
    cout<<X[i]<<" "<<EXP[i]<<" "<<EXPM1[i]<<" "<<EXPP1[i]<<" "<<EXPM2[i]<<" "<<EXPP2[i]<<endl;
    xmass =X[i];
    i++ ;
    if (xmass==maxmass) break;
  }

  expFile.close();

  ifstream obsFile;
  obsFile.open(obsfile.c_str());

  i=0;
  xmass=0;
  while( !obsFile.eof() ){
    obsFile >> X[i] >> OBS[i];
    xmass =X[i];
    cout << i <<" " << xmass << " " <<OBS[i] <<endl;
    i++ ;
    if (xmass==maxmass) break;
  }

  obsFile.close();

  int Np=i-1;

  ifstream theoryFile;
  theoryFile.open(theoryfile.c_str());

  i=0;
  xmass=0.;
  while(!theoryFile.eof() ){
    theoryFile >> X[i] >> CS[i];
    CS[i]=CS[i]*1.2;
    CS2[i]=CS[i]*4/3;

    xmass =X[i];
    i++ ;
    //cout << i <<" " << xmass<<endl;
    if (xmass==maxmass) break;

  }

  int Ndata=i;
  cout << "Ndata " << Ndata << endl;

  cout<< "mass"<<" theory "<<" observed"<<" +2 sigma"<<" +1 sigma"<<" expected"<<"  -1 sigma"<<" -2 sigma"<<endl;

  for (int i=0; i<Ndata-1;i++) {
    cout<< X[i]<< "  " << CS[i]<< "  "<< " " << OBS[i]<< "   "<<"  "<<EXPP2[i]<<"   "<<EXPP1[i]<<"  "<<EXP[i]<<"  "<<EXPM1[i]<<"  "<<EXPM2[i]<<endl;
  }

  TGraph *xsec = new TGraph(Ndata,X,CS);
  xsec->SetName("xsec");
  TGraph *xsec2 = new TGraph(Ndata,X,CS2);
  xsec2->SetName("xsec2");
  TGraph *lobs = new TGraph(Ndata,X,OBS);
  lobs->SetName("lobs");
  TGraph *lexp = new TGraph(Ndata,X,EXP);
  lexp->SetName("lexp");

  xsec->SetTitle("");
  xsec->SetLineColor(kRed);

  xsec2->SetTitle("");
  xsec2->SetLineColor(kRed);

  lobs->SetLineColor(kBlack);
  lexp->SetLineColor(kBlack);

  xsec->SetLineWidth(3);
  xsec2->SetLineWidth(3);
  xsec2->SetLineStyle(5);

  lobs->SetLineWidth(2);
  lexp->SetLineWidth(2);

  lexp->SetLineStyle(2);

  if (maxmass==2500) {
    float xbins[15]={750,850,950,1050,1150,1250,1350,1450,1550,1650,1700,2000,2200,2400,2600};
  };
  if (maxmass==2300) {
    float xbins[15]={750,850,950,1050,1150,1250,1350,1450,1550,1650,1700,2000,2200,2400};
  };
  if (maxmass==2100) {
    float xbins[15]={750,850,950,1050,1150,1250,1350,1450,1550,1650,1700,2000,2200};
  };
  if (maxmass==1900) { 
    float xbins[12]={750,850,950,1050,1150,1250,1350,1450,1550,1650,1750,2050};
  };
  if (maxmass==3000) { 
    float xbins[21]={1050,1150,1250,1350,1450,1550,1650,1750,1850,1950,2050,2150,2250,2350,2450,2550,2650,2750,2850,2950,3050};
  }

  TH1D *hr = new TH1D("hr","",Ndata,xbins);

  //  for(int i=0;i<Ndata;i++) {
  //    float center = (EXPP2[i]+EXPM2[i])/2.0;
  //    hrx->SetBinContent(i+1,center);
  //    hrx->SetBinError(i+1,TMath::Abs(EXPM2[i]-center));
  //    hrerr->SetBinContent(i+1,TMath::Abs(EXPM2[i]-center)/center);
  //  }

  for(int i=0;i<Ndata;i++) {
    float center = (EXPP2[i]+EXPM2[i])/2.0;
    float error = TMath::Abs(EXPM2[i]-center);
    hr->SetBinContent(i+1,center);
    hr->SetBinError(i+1,error);
  }

  TH1D *hr2 = new TH1D("hr2","",Ndata,xbins);

  //  for(int i=0;i<Ndata;i++) {
  //    float center = (EXPP1[i]+EXPM1[i])/2.0;
  //    hrx->SetBinContent(i+1,center);
  //    hrx->SetBinError(i+1,TMath::Abs(EXPM1[i]-center));
  //    hrerr->SetBinContent(i+1,TMath::Abs(EXPM1[i]-center)/center);
  //  }

  for(int i=0;i<Ndata;i++) {
    float center = (EXPP1[i]+EXPM1[i])/2.0;
    float error = TMath::Abs(EXPM1[i]-center);
    hr2->SetBinContent(i+1,center);
    hr2->SetBinError(i+1,error);
  }

  TCanvas *c1 = new TCanvas("c1","PRL Fig. 3",1000,800);  
  c1->cd();
  gStyle->SetOptStat(0);
  c1->Range(0,0,1,1);
  c1->SetTopMargin(0.06);

  hr->SetFillColor(kYellow);
  hr->SetLineColor(kBlack);
  hr->SetMarkerColor(0);
  hr->SetMarkerSize(0);
  hr->SetLineWidth(0);
  if (limitFILE.find(sRH)<30) {
    hr->GetYaxis()->SetTitle("#sigma(pp#rightarrowW'_{R}) #times B(W'_{R}#rightarrowtb#rightarrowl#nubb) [pb]");
    hr->GetXaxis()->SetTitle("W'_{R} Mass [GeV]");
  }
  else {
    hr->GetYaxis()->SetTitle("#sigma(pp#rightarrowW'_{R}) #times B(W'_{R}#rightarrowtb#rightarrowl#nubb) [pb]");
    hr->GetXaxis()->SetTitle("W' Mass [GeV]");
  }
  hr->SetMaximum(20);
  hr->SetMinimum(0.004);
  if (limitFILE.find(sRH)<30) {    
    hr->SetMaximum(6.0);
    hr->SetMinimum(0.004);
  }
  hr->SetMaximum(16.0);
  hr->SetMinimum(0.0004);

  hr->Draw("E3");

  hr2->SetFillColor(kGreen);
  hr2->SetLineColor(kBlack); 
  hr2->SetLineWidth(0); 
  hr2->SetMarkerColor(0);   
  hr2->SetMarkerSize(0);   
  hr2->GetXaxis()->SetLabelOffset(0.0001);
  hr2->Draw("E3same");

  lexp->Draw("L");
  xsec->Draw("C");
  xsec2->Draw("C");

  c1->SetLogy();
    
  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(31);// # align right
  latex2->DrawLatex(0.87, 0.95, "CMS Preliminary, 12.2 fb^{-1} at #sqrt{s} = 8 TeV");

  TLatex *latex3 = new TLatex();
  latex3->SetNDC();
  latex3->SetTextSize(0.05);
  latex3->SetTextAlign(22); // # align right
  latex3->DrawLatex(0.37, 0.20, leptonType);

  TLatex *latex4 = new TLatex();
  latex4->SetNDC();
  latex4->SetTextSize(0.05);
  latex4->SetTextAlign(22); // # align right
  latex4->DrawLatex(0.38, 0.27, "Invariant Mass Analysis");
    
  TLegend *leg = new TLegend(.51,.56,.9,.91,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(10);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  TLegendEntry *entry=leg->AddEntry("xsec","Theory M_{#nu_{R}}<< M_{W'_{R}}","l");
  TLegendEntry *entry=leg->AddEntry("xsec2","Theory M_{#nu_{R}}> M_{W'_{R}}","l");
  entry=leg->AddEntry("lobs","95% CL observed","lp");
  entry=leg->AddEntry("lexp","95% CL expected","l");
  entry=leg->AddEntry(hr2,"#pm1#sigma expected","f");
  entry=leg->AddEntry(hr,"#pm2#sigma expected","f");
  leg->Draw();

  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  char pname[100];
  sprintf(pname,"./limit_plots/%s_bayesian_12pt2.png",limitFILE.c_str());
  c1->Print(pname);
  sprintf(pname,"./limit_plots/%s_bayesian_12pt2.pdf",limitFILE.c_str());
  c1->Print(pname);
  sprintf(pname,"./limit_plots/%s_bayesian_12pt2.root",limitFILE.c_str());
  c1->Print(pname);

}
