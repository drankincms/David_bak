{
  //#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
  //#include<stdio.h>

using namespace std;


//"inner" array go through ar values. each inner array corresponds to al value. in  exp_1000[i][j] i goes through al values, j ar values.
// 21 mass values
 Double_t obslim[11][11][21];
 Double_t explim[11][11][21];
 Double_t theory[11][11];


//open file to save histograms
TFile *file=new TFile("plots_observed.root","recreate");
gStyle->SetTextFont(22);

//make a 3d array (ar,al,mass) 
Double_t x[11][11][23],t[11][11][23];
 Double_t m[21]={800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,2800,2900};
Double_t ar[11]={0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
Double_t al[11]={0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
Double_t mlim[11][11];
Double_t mlim_al[21][11];

//Double_t theoryR[23]= {2.3022,1.3818,0.85538,0.54325,0.35203,0.23219,0.15547,0.10518,0.072012,0.049683,0.034576,0.024249,0.017124,0.012176,0.0087191,0.0062918,0.0045757,0.0033568,0.0024870,0.0018624,0.0014102,0.0010818,0.00084115};

Double_t theoryR[21]= {2.3022,1.3818,0.85538,0.54325,0.35203,0.23219,0.15547,0.10518,0.072012,0.049683,0.034576,0.024249,0.017124,0.012176,0.0087191,0.0062918,0.0045757,0.0033568,0.0018624,0.0014102,0.0010818};

// Double_t theoryL[23] = {3.1080,2.2731,1.8087,1.547,1.3870,1.2945,1.2390,1.2061,1.1869,1.1761,1.1705,1.1678,1.1673,1.1680,1.1692,1.1711,1.1727,1.1746,'2600':1.1763,1.1780,1.1797,1.1810,1.1825};

Double_t theoryL[21] = {3.1080,2.2731,1.8087,1.547,1.3870,1.2945,1.2390,1.2061,1.1869,1.1761,1.1705,1.1678,1.1673,1.1680,1.1692,1.1711,1.1727,1.1746,1.1780,1.1797,1.1810};

// Double_t theoryM[23]= {5.4166,3.6684,2.6815,2.1031,1.7539,1.5389,1.4043,1.3194,1.2650,1.2305,1.2090,1.1954,1.1872,1.1824,1.1798,1.1787,1.1784,1.1791,'2600':1.1792,1.1803,1.1813,1.1825,1.1835};

Double_t theoryM[21]= {5.4166,3.6684,2.6815,2.1031,1.7539,1.5389,1.4043,1.3194,1.2650,1.2305,1.2090,1.1954,1.1872,1.1824,1.1798,1.1787,1.1784,1.1791,1.1803,1.1813,1.1825};


Double_t theorySM=5.55; // in pb




	
bool set=false;
//Wprime_Histos_combined-4-1
const char base[] = "run_analysis/Wprime_Histos_combined-";
const char mid[] = "-";
char filename [80];
 char line[80];
string header;
float mx[25],y[25],y0lo[25],y0hi[25],y1lo[25],y1hi[25];
 int num=0;
 int numal=0;
 int numar=0;

 for (int il=0; il<11; il++){
   for (int ir=0; ir<11; ir++){
     //if (ir==0 && il==0){ continue;}
     //       cout << il << " "<< ir << endl;

     //       cout << il << " "<< ir << endl;
     numal=il;
     numar=ir;
     sprintf(filename, "%s%d%s%d%s", base, numal,mid,numar,"/bayesian_limits_observed.txt");
     cout << filename << endl;
     num++;
     //       cout << " analysis num " << num << endl;
     ifstream inFile;
     inFile.open(filename);
     for (int i=0;i<3;i++){
       inFile >> line;
       //cout << line << endl;
     }
     int i=0;
     while(!inFile.eof() ){
       inFile >> mx[i] >> y[i];
       cout<<mx[i]<<" "<<y[i] <<endl;
       if (i<21) {
         explim[il][ir][i]=y[i];
         x[il][ir][i]=y[i];
       }
       i++;
     }
     int imax=i-1;
     for (int ip=0;ip<imax;ip++){
       //      cout<<mx[ip]<<" "<<y[ip]<<" "<<y0lo[ip]<<" "<<y0hi[ip]<<" "<<y1lo[ip]<<" "<<y1hi[ip]<<endl;
     }
   }

 }


 // cout << il << " "<< ir << endl;
 /// this is the theory cross section for each al/ar combination -
 for (int k=0; k<11; k++) {
   for (int j=0; j<11; j++) {
     double aal=0.1*k;
     double aar=0.1*j;
     //     cout << " al/ar " << aal << " " << aar << k << " " << j <<endl;
     double rt_wgt=(-aal*aal + aal*aal*aal*aal + aar*aar*aar*aar- aal*aal*aar*aar);
       double lt_wgt=(aal*aal - aal*aal*aar*aar);
       double st_wgt=(1-aal*aal);
       double mt_wgt=(aal*aal*aar*aar);
       for (int im=0; im<21; im++){
	 t[k][j][im]=(1.2*(rt_wgt*theoryR[im]+lt_wgt*theoryL[im]+mt_wgt*theoryM[im])+st_wgt*theorySM);
	 //x[k][j][im]=x[k][j][im];
	 x[k][j][im]=x[k][j][im]*t[k][j][im];
	 cout << "al/ar/ mass "<< " " << aal << " " << aar << " " << m[im] << " " << x[k][j][im]<< " " << t[k][j][im]<< endl; 

       }
   }
 }

    


 //========================.
 // below this is the D0 code. there we had 9 masses... now we have 21 or 23  masses for CMS 
 //

 for (int i=0;i<11;i++){
   for (int j=0;j<11;j++){
     mlim[i][j]=-1;
     //find intersection between observed limit and predicted cross section
     for (int k=0;k<20;k++){
       Double_t xxx=(x[i][j][k+1]-t[i][j][k+1])*(x[i][j][k]-t[i][j][k]);
       //cout<<x[i][j][k+1]<<" "<<t[i][j][k+1]<<" "<<xxx<<endl;
       if (xxx<0){
	 Double_t f=(x[i][j][k]-t[i][j][k])/((x[i][j][k]-t[i][j][k])-(x[i][j][k+1]-t[i][j][k+1]));
	 if ( mlim[i][j] = -1) {
	 mlim[i][j]=m[k]+(m[k+1]-m[k])*f;
	 cout<< "al/ar/explo/exphi/thlo/explo/mlim " << i << " " << j << " " << x[i][j][k] << " " << x[i][j][k+1] << " " << t[i][j][k] << " " << t[i][j][k+1]<< " " << mlim[i][j]<<endl;
	 set=true;
	 }
	 else {
	   set=false;
	 }
	 //cout<<mlim[i][j]<<endl;    
       }//if
     }//k
     if (!set) mlim[i][j]=799; set=false;
   }//j
 }// i
    
 //make a contour plot from the original grid ar vs al
 TH2 *obs2 = new TH2F("obs2","ar vs al", 11,-0.05,1.05, 11,-0.05, 1.05);
 for (int i=0;i< 11;i++){
   for (int j=0;j< 11;j++){
     //cout << mlim[i][j] << "; " ;
     obs2->SetBinContent(i+1,j+1,mlim[i][j]);
   }
   //cout <<endl;
 }
    
 Double_t ax[24],at[24];

 //interpolate linearly between grid points to make a finer grid
 Double_t marray[101][101];
 for (int j=0;j<101;j++){
   float s=j;
   s=s/100;
   int js=int(s*10);
   //cout<<js<<endl;
   for (int i=0;i<101;i++){
     float r=i;
     r=r/100;
     int ir=int(r*10);
     //cout << ir<<endl;
     for (int k=0;k<21;k++){
       if (ir<10){
	 float ax1=x[ir][js][k]+(10*r-ir)*(x[ir+1][js][k]-x[ir][js][k]);
	 float at1=t[ir][js][k]+(10*r-ir)*(t[ir+1][js][k]-t[ir][js][k]);
       }else{
	 float ax1=x[10][js][k];
	 float at1=t[10][js][k];
       }
       ax[k]=ax1;
       at[k]=at1;
       if (js<10){
	 if (ir<10){
	   float ax2=x[ir][js+1][k]+(10*r-ir)*(x[ir+1][js+1][k]-x[ir][js+1][k]);
	   float at2=t[ir][js+1][k]+(10*r-ir)*(t[ir+1][js+1][k]-t[ir][js+1][k]);
	 }else{
	   float ax2=x[10][js+1][k];
	   float at2=t[10][js+1][k];
	 }
       }
       ax[k]=ax1+(10*s-js)*(ax2-ax1);
       at[k]=at1+(10*s-js)*(at2-at1);
       //cout << ax[k] <<";"<< at[k] << endl;
     }

     //find the mass limit
     //cout<<ax[0]<<" "<<at[0]<<endl;
     for (int k=0;k<20;k++){
       marray[i][j]=799;
       Double_t xxx=(ax[k+1]-at[k+1])*(ax[k]-at[k]);
       //cout<<ax[k+1]<<" "<<at[k+1]<<" "<<xxx<<endl;
       if (xxx<0){

	 Double_t f=(ax[k]-at[k])/((ax[k]-at[k])-(ax[k+1]-at[k+1]));
	 marray[i][j]=m[k]+(m[k+1]-m[k])*f;
	 //cout<<i<<" "<<j<<" "<<marray[i][j]<<endl;
	 break;
       }
     }
   }
 }

 //now we fill a histogram with the finer grid

 TH2 *obs1 = new TH2F("obs1","ar vs al", 101,-0.005,1.005, 101,-0.005, 1.005);
 for (int i=0;i<101;i++){
   for (int j=0;j<101;j++){
     //cout << marray[i][j] << ";" ;
     //mass limit cannot be below 800 because we have no data below 800
     if (marray[i][j]<799) marray[i][j]=799;
     obs1->SetBinContent(i+1,j+1,marray[i][j]);
   }
   //cout <<endl;
 }
   
 TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
color=(TColor*)(gROOT->GetListOfColors()->At(100));
color->SetRGB(1,0,0);
color=(TColor*)(gROOT->GetListOfColors()->At(99));
color->SetRGB(1,0.1,0.1);
color=(TColor*)(gROOT->GetListOfColors()->At(98));
color->SetRGB(1,0.2,0.2);
color=(TColor*)(gROOT->GetListOfColors()->At(97));
color->SetRGB(1,0.3,0.3);
color=(TColor*)(gROOT->GetListOfColors()->At(96));
color->SetRGB(1,0.4,0.4);
color=(TColor*)(gROOT->GetListOfColors()->At(95));
color->SetRGB(1,0.5,0.5);
color=(TColor*)(gROOT->GetListOfColors()->At(94));
color->SetRGB(1,0.6,0.6);
color=(TColor*)(gROOT->GetListOfColors()->At(93));
color->SetRGB(1,0.7,0.7);
color=(TColor*)(gROOT->GetListOfColors()->At(92));
color->SetRGB(1,0.8,0.8);
color=(TColor*)(gROOT->GetListOfColors()->At(91));
color->SetRGB(1,0.9,0.9);
color=(TColor*)(gROOT->GetListOfColors()->At(90));
color->SetRGB(1,1,1);
Int_t color_indices[11]={90,91,92,93,94,95,96,97,98,99,100};
gStyle->SetPalette(11,color_indices);

 c1->SetFillColor(kWhite); //add                                                                                           
 c1->SetBottomMargin(0.15); //change 2 0.15
 c1->SetRightMargin(0.15);
 c1->SetLeftMargin(0.17);
 obs1->SetContour(20,m);

 obs1->GetZaxis()->SetLabelFont(22);

// obs1->Draw("cont4z");
 obs1->Draw("cont3zcol");
 obs1->SetLineWidth(3);
 c1->SetBottomMargin(0.15); //change 2 0.15
 c1->SetRightMargin(0.15);
 c1->SetLeftMargin(0.17);
 obs1->SetContour(20,m);
 obs1->Draw("cont2same");  
 obs1->GetYaxis()->SetTitleFont(22);
 obs1->GetXaxis()->SetTitleFont(22);
 obs1->GetYaxis()->SetTitleSize(0.07);
 obs1->GetXaxis()->SetTitleSize(0.07);
 obs1->GetXaxis()->SetLabelFont(22); //add times new roman 'bold'
 obs1->GetYaxis()->SetLabelFont(22); //add times new roman 'bold' 
 obs1->GetXaxis()->SetTitle("a^{L}");
 obs1->GetYaxis()->SetTitle("a^{R}");
 //obs1->GetXaxis()->SetTitleOffset(-0.5); //add


TLatex *t1 = new TLatex();
float x1=0.2;
float y1=.85;

t1->SetNDC();
t1->SetTextFont(22);
t1->SetTextColor(1);
t1->SetTextAlign(11);
t1->SetTextSize(0.06);
t1->SetTextFont(22);

TLatex *t1temp=0;
//t1temp=t1->DrawLatex(x1,y1,"(b)   D\351 2.3 fb^{  -1}");
//t1temp=t1->DrawLatex(x1,y1-0.09,"Limits on M(W')");
c1->Update();
c1->Print("contour_arVSal-combined_observed.pdf");
c1->Print("contour_arVSal-combined_observed.png");

// TCanvas *c2 = new TCanvas("c2","c2",600,0,600,600);

// obs1->Draw("text");
 obs1->Write();
 obs2->Write();
    
}//file
  
