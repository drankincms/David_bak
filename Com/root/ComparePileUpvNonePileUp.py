#!/usr/bin/env python

import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

from tdrStyle import *
setTDRStyle()


def cmsPrel(intLumi, x=0.922771, y=0.86756):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)

    latex.SetTextAlign(31) # align right
    latex.DrawLatex(x, y,"#sqrt{s} = 7 TeV")
    if (intLumi > 0.):
        latex.SetTextAlign(31) # align right
        #latex.DrawLatex(x+0.02,y-0.08,"#int #font[12]{L}dt = "+str(intLumi)+"nb^{-1}")
        latex.DrawLatex(x+0.02, y-0.08,"#int #font[12]{L}dt = "+'%.1f'%intLumi+"pb^{-1}")
    
    latex.SetTextAlign(5) # align left
    #latex.DrawLatex(0.160643,0.96131,"#font[22]{CMS preliminary 2010}")


Data = TFile("EWK_met_W/Data/All.root")
MCWmunu = TFile("EWK_met_W/Wmunu_pythia38/All.root")
MCWmunu_PileUp = TFile("EWK_met_W/Wmunu_pythia38_PileUp/All.root")

TData = Data.Get("mySUperTree"); 
TMCWmunu = MCWmunu.Get("mySUperTree"); 
TMCWmunu_PileUp = MCWmunu_PileUp.Get("mySUperTree"); 

cut = ""
cut1 = "muon_pt > 25 && PF_met_pt > 25 && PF_Mt > 50"
cut2 = "muon_pt > 25 && PF_met_pt > 25 && PF_Mt > 50 && nHitsPV == 1"
cut3 = "muon_pt > 25 && PF_met_pt > 25 && PF_Mt > 50 && nHitsPV > 1"


##### Pick you bins
bins = 60
low = 0
high = 120

##### Which cut version do you want
cuttype = cut2

###### What variable do you want to run over

#var = "nHitsPV"
#varlabel = "Number of Primary Verticies"
var = "PF_met_pt"
varlabel = "#slash{E}_{T}"
#var = "PF_Mt"
#varlabel = "M_{T}"
#var = "PF_recoil"
#varlabel = "u_{T}"

###### Are you requiring 1 or All PV
Pvlabel = "Reqruie 1 PV"
save = var + "_1PV"
#Pvlabel = "All PV"
#save = var + "_AllPV"


Data_PF_met_pt = TH1D("Data_PF_met_pt", "#slash{E}_{T}", bins, low, high);
TData.Project("Data_PF_met_pt", var, cuttype);

Data_PF_met_pt_ratio = TH1D("Data_PF_met_pt_ratio", "#slash{E}_{T}", bins, low, high);
TData.Project("Data_PF_met_pt_ratio", var, cuttype);

Wmunu_PF_met_pt = TH1D("Wmunu_PF_met_pt", "#slash{E}_{T}", bins, low, high);
TMCWmunu.Project("Wmunu_PF_met_pt", var, cuttype);

Wmunu_PF_met_pt_ratio = TH1D("Wmunu_PF_met_pt_ratio", "#slash{E}_{T}", bins, low, high);
TMCWmunu.Project("Wmunu_PF_met_pt_ratio", var, cuttype);

Wmunu_PF_met_pt_PileUp = TH1D("Wmunu_PF_met_ptPileUp", "#slash{E}_{T}", bins, low, high);
TMCWmunu_PileUp.Project("Wmunu_PF_met_ptPileUp",var, cuttype);

Wmunu_PF_met_pt.Scale(Data_PF_met_pt.Integral() / Wmunu_PF_met_pt.Integral() ); 
Wmunu_PF_met_pt_ratio.Scale(Data_PF_met_pt.Integral() / Wmunu_PF_met_pt_ratio.Integral() ); 
Wmunu_PF_met_pt_PileUp.Scale(Data_PF_met_pt.Integral() / Wmunu_PF_met_pt_PileUp.Integral() ); 

Wmunu_PF_met_pt_ratio.Divide(Wmunu_PF_met_pt_PileUp)
Data_PF_met_pt_ratio.Divide(Wmunu_PF_met_pt_PileUp)

c1 = TCanvas("c1","Full PF Cuts for Data vs WJets MC", 1000, 800);
c1.Divide(2,2) 
c1.cd(1) 
Wmunu_PF_met_pt.Draw()
Wmunu_PF_met_pt.GetXaxis().SetTitle("No Pile Up " + varlabel  +", "+ Pvlabel);

leg = TLegend(0.5548779,0.622797,0.932312,.927762);
leg.AddEntry( Wmunu_PF_met_pt, "W #rightarrow #mu #nu No Pile Up","l");
leg.SetShadowColor(0);
leg.SetFillColor(0);
leg.SetLineColor(0);
leg.Draw();

c1.cd(2) 
Wmunu_PF_met_pt_PileUp.Draw()
Wmunu_PF_met_pt_PileUp.GetXaxis().SetTitle("With Pile Up " + varlabel  +", "+ Pvlabel);

leg2 = TLegend(0.5548779,0.622797,0.932312,.927762);
leg2.AddEntry( Wmunu_PF_met_pt_PileUp, "W #rightarrow #mu #nu With Pile Up","l");
leg2.SetShadowColor(0);
leg2.SetFillColor(0);
leg2.SetLineColor(0);
leg2.Draw();


c1.cd(3) 
Data_PF_met_pt.SetMaximum(Data_PF_met_pt.GetMaximum() + 0.2*Data_PF_met_pt.GetMaximum())
Data_PF_met_pt.Draw("E1")
Wmunu_PF_met_pt.Draw("SAME")
Wmunu_PF_met_pt_PileUp.SetLineColor(2)
Wmunu_PF_met_pt_PileUp.Draw("SAME")
Data_PF_met_pt.GetXaxis().SetTitle("Compare " + varlabel  +", "+ Pvlabel);

leg1 = TLegend(0.5548779,0.622797,0.932312,.927762);
leg1.AddEntry( Data_PF_met_pt, "Data", "lp");
leg1.AddEntry( Wmunu_PF_met_pt, "W #rightarrow #mu #nu No Pile Up","l");
leg1.AddEntry( Wmunu_PF_met_pt_PileUp, "W #rightarrow #mu #nu With Pile Up","l");
leg1.SetShadowColor(0);
leg1.SetFillColor(0);
leg1.SetLineColor(0);
leg1.Draw();

cmsPrel(8.57, 0.822771, 0.46756)


c1.cd(4)
Wmunu_PF_met_pt_ratio.Draw("E1")
Wmunu_PF_met_pt_ratio.GetXaxis().SetTitle(varlabel  +", "+ Pvlabel);
Wmunu_PF_met_pt_ratio.GetYaxis().SetTitle("No Pile Up / Pile Up");



c1.SaveAs("ComparePileUp_"+save+".gif")





c2 = TCanvas("c2","Full PF", 1000, 800);
c2.Divide(1,2) 
c2.cd(1) 

Data_PF_met_pt.SetMaximum(Data_PF_met_pt.GetMaximum() + 0.2*Data_PF_met_pt.GetMaximum())
Data_PF_met_pt.Draw("E1")
Wmunu_PF_met_pt_PileUp.SetLineColor(2)
Wmunu_PF_met_pt_PileUp.Draw("SAME")
Data_PF_met_pt.GetXaxis().SetTitle( varlabel );
    
leg11 = TLegend(0.5548779,0.622797,0.932312,.927762);
leg11.AddEntry( Data_PF_met_pt, "Data", "lp");
leg11.AddEntry( Wmunu_PF_met_pt_PileUp, "W #rightarrow #mu #nu With Pile Up","l");
leg11.SetShadowColor(0);
leg11.SetFillColor(0);
leg11.SetLineColor(0);
leg11.Draw();

cmsPrel(8.57, 0.822771, 0.46756)

c2.cd(2)
Data_PF_met_pt_ratio.Draw("E1")
Data_PF_met_pt_ratio.GetXaxis().SetTitle(varlabel);
Data_PF_met_pt_ratio.GetYaxis().SetTitle("Data / Pile Up");
        
c2.SaveAs("ComparePileUp_DatatoPileUp_"+save+".gif")







