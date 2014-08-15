import ROOT, sys, os, string, re
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2D, TH2F, TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack, TGraph, TGraphErrors, TColor, TMath
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad

from array import array

import math
from math import *

from tdrStyle import *
setTDRStyle()

from LoadData import *

global ksfile
ksfile = open('Wprime_Plots/VARKSCHI2val.txt', 'w')

def plot_forPAS(channel, var, bin, low, high, ylabel, xlabel, save, nBtags = -1, setLog = False, finalcuts = False):

    if (channel == 'electron'):
        List = ['Data_el','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph','Wprime1700Right','Wprime1900Right','Wprime2100Right']
        #List = ['Data_el','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']
    if (channel == 'muon'):
        List = ['Data_mu','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph','Wprime1700Right','Wprime1900Right','Wprime2100Right']
        #List = ['Data_mu','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']
      

    if (channel == 'electron'):      
        #cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 32 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 32 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 


    if (channel == 'muon'):      
        #cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20'
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20'


    if finalcuts: cut+= '&& BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 75  && Jet1Jet2_Pt_LjetsTopoCalcNew > 100'  
 
    njets = ""

    if nBtags == -1:
        cutbtag = ''
        save = save     
    if nBtags == 0:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) == 0 )'
        save = save + '_0bTags'
        njets = "N_{b tags} = 0"
    if nBtags == 1:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) == 1)'
        save = save + '_1bTags'
        njets = "N_{b tags} = 1"
    if nBtags == 2:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) >= 1 ) '
        save = save + '_GE1bTags'
        njets = "N_{b tags} #geq 1"
    if nBtags == 3:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc==1 && (jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) || (jet_1_tag_WprimeCalc==1 && (jet_0_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) ) '
        save = save + '_GE2bTags'
        njets = "N_{b tags} #geq 2"


    #cutwjets = ' && nSelJets_CommonCalc <3 && Ht_LjetsTopoCalcNew < 300 '
    #cutttbar = ' && nSelJets_CommonCalc>4 '
 
    cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
    cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
    cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light
    #cutwbb = ' ' # Wb(b)
    #cutwcc = ' ' # Wc(c)
    #cutwjj = ' ' # W+light
                                                                                                                   
    #SFWjmu = 0.85*1.0
    #SFWcmu = 0.92*1.66
    #SFWbmu = 0.92*1.21
    #SFWjmu = 1.0
    #SFWcmu = 1.0
    #SFWbmu = 1.0
    #SFWjmu = 0.9*1.0
    #SFWcmu = 1.15*1.66
    #SFWbmu = 1.15*1.21
    #SFWjmu = 1.13*0.85
    #SFWcmu = 1.24*0.92*1.66
    #SFWbmu = 1.24*0.92*1.21
    #SFWjmu = 1.03*0.85
    #SFWcmu = 1.27*0.92*1.66
    #SFWbmu = 1.27*0.92*1.21
    #SFWjmu = 1.05*0.85       ## myHF
    #SFWcmu = 1.25*0.92*1.66  ## myHF
    #SFWbmu = 1.25*0.92*1.21  ## myHF
    SFWjmu = 1.08*0.85      ## myHF120
    SFWcmu = 1.06*0.92*1.66  ## myHF120
    SFWbmu = 1.06*0.92*1.21  ## myHF120

    save = save + '_myHF_jet0pt120'

    WjjHist = TH1D('WjjHist', 'WjjHist', bin,low,high)
    WccHist = TH1D('WccHist', 'WccHist', bin,low,high)
    WbbHist = TH1D('WbbHist', 'WbbHist', bin,low,high)

    WjjHistPre = TH1D('WjjHistPre', 'WjjHistPre', bin,low,high)
    WccHistPre = TH1D('WccHistPre', 'WccHistPre', bin,low,high)
    WbbHistPre = TH1D('WbbHistPre', 'WbbHistPre', bin,low,high)

    Variables = {}
    VariablesPre= {}
    efficiency = {}
    
    if (channel=='electron'):
        weight = 'weight_PU_ABC_PileUpCalc*weight_ElectronEff_WprimeCalc'
    if (channel=='muon'):
        weight = 'weight_PU_ABC_PileUpCalc*weight_MuonEff_WprimeCalc'

    for Type in List:
        Variables[Type] = TH1D(Type+var+channel, Type+var+channel, bin, low, high)
        VariablesPre[Type] = TH1D(Type+var+channel+'Pre', Type+var+channel+'Pre', bin, low, high)

        histName = Type+var+channel
        histNamePre = Type+var+channel+'Pre'
     
        if Type == 'WJets':
            #print 'Filling ',Type
            Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutbtag + cutwjj + ")", 'goff')
            Trees[Type].Draw(var + " >> " + "WbbHist","("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutbtag + cutwbb + ")", 'goff')
            Trees[Type].Draw(var + " >> " + "WccHist","("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutbtag + cutwcc + ")", 'goff')
            #print 'Raw Wjj ',Variables['WJets'].Integral()
            Variables[Type].Add(WbbHist)
            #print 'Raw Wjj + Wbb ',Variables['WJets'].Integral()
            Variables[Type].Add(WccHist) 
            #print 'Raw Wjj + Wbb + Wcc',Variables['WJets'].Integral()
            Trees[Type].Draw(var + " >> " + histNamePre, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
            Trees[Type].Draw(var + " >> " + "WbbHistPre","("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
            Trees[Type].Draw(var + " >> " + "WccHistPre","("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
            #print 'Pretag Raw Wjj ',VariablesPre['WJets'].Integral()
            VariablesPre[Type].Add(WbbHistPre)
            #print 'Pretag Raw Wjj + Wbb ',VariablesPre['WJets'].Integral()
            VariablesPre[Type].Add(WccHistPre) 
            #print 'Pretag Raw Wjj + Wbb + Wcc ',VariablesPre['WJets'].Integral()

        elif Type.startswith('Data'):
            #print 'Filling ',Type
            Trees[Type].Draw(var + " >> " + histName, "(" + cut + cutbtag + ")", 'goff')
        else:
            #print 'Filling ',Type
            Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + cut + cutbtag + ")", 'goff')
            Trees[Type].Draw(var + " >> " + histNamePre, "("+weight+")*(" + cut + ")", 'goff')
            
        if (not Type.startswith('Data')):

            SF = 1.0
            if (channel == 'electron'):
                lumi = lumi_el
                SF = 0.977
            if (channel == 'muon'):
                lumi = lumi_mu
                #SF = 0.966

            if (Type.startswith('Wprime')): SF *= 20
            #print 'EVENTS Before Scaling FOR ',Type,' = ',Variables[Type].Integral()
            #print 'Pre Events before scaling for ',Type,' = ',VariablesPre[Type].Integral()
            #print str(SF),' ',str(lumi),' ',str(xsec_norm[Type]),' ',str(Nevents[Type])

            if Variables[Type].Integral() != 0:
                #print Type,' Lumi scaling: ',str(SF*lumi*xsec[Type]/Nevents[Type])
                Variables[Type].Scale ( SF*lumi*xsec[Type]/Nevents[Type] ) 
                VariablesPre[Type].Scale ( SF*lumi*xsec[Type]/Nevents[Type] ) 
                if ( (not Type.startswith('T')) and (not Type.startswith('Wprime')) ):
                    dummy = 1.0
                    #print Type,' Pretag scaling: ',str(Variables[Type].Integral() / VariablesPre[Type].Integral())
                    #VariablesPre[Type].Scale (  Variables[Type].Integral() / VariablesPre[Type].Integral() )
                    #Variables[Type] = VariablesPre[Type].Clone()
                efficiency[Type] = Variables[Type].Integral()/Nevents[Type]
            else:
                efficiency[Type] = 0

            #Variables[Type].SetBinContent(bin,Variables[Type].GetBinContent(bin)+Variables[Type].GetBinContent(bin+1) ) 
            #Variables[Type].SetBinContent(bin+1,0) 
            #VariablesPre[Type].SetBinContent(bin,VariablesPre[Type].GetBinContent(bin)+VariablesPre[Type].GetBinContent(bin+1) ) 
            #VariablesPre[Type].SetBinContent(bin+1,0) 

            #print 'SCALED EVENTS FOR ',Type,'  = ',Variables[Type].Integral()
    
    Variables['TTbar_Madgraph'].SetFillColor(ROOT.kRed-7)
    Variables['T_s'].SetFillColor(ROOT.kRed-7)
    Variables['Tbar_s'].SetFillColor(ROOT.kRed-7)
    Variables['T_t'].SetFillColor(ROOT.kRed-7)
    Variables['Tbar_t'].SetFillColor(ROOT.kRed-7)
    Variables['T_tW'].SetFillColor(ROOT.kRed-7)
    Variables['Tbar_tW'].SetFillColor(ROOT.kRed-7)

    Variables['T_s'].SetLineColor(ROOT.kRed-7)
    Variables['Tbar_s'].SetLineColor(ROOT.kRed-7)
    Variables['T_t'].SetLineColor(ROOT.kRed-7)
    Variables['Tbar_t'].SetLineColor(ROOT.kRed-7)
    Variables['T_tW'].SetLineColor(ROOT.kRed-7)
    Variables['Tbar_tW'].SetLineColor(ROOT.kRed-7)

    Variables['WJets'].SetFillColor(ROOT.kGreen-3)
    Variables['ZJets_M50'].SetFillColor(ROOT.kGreen-3)
    Variables['WW'].SetFillColor(ROOT.kGreen-3)
    Variables['ZJets_M50'].SetLineColor(ROOT.kGreen-3)
    Variables['WW'].SetLineColor(ROOT.kGreen-3)
    
    Variables['TTbar_Madgraph'].SetLineWidth(2)
    Variables['WJets'].SetLineWidth(2)

    #Variables['Wprime1700Right'].SetLineColor(1)
    #Variables['Wprime1700Right'].SetLineWidth(2)
    #Variables['Wprime1700Right'].SetLineStyle(4)

    #Variables['Wprime1900Right'].SetLineColor(1)
    #Variables['Wprime1900Right'].SetLineWidth(2)
    #Variables['Wprime1900Right'].SetLineStyle(6)

    #Variables['Wprime2100Right'].SetLineColor(1)
    #Variables['Wprime2100Right'].SetLineWidth(2)
    #Variables['Wprime2100Right'].SetLineStyle(8)

    stack = THStack('a', 'a')
    added = TH1D('a', 'a',bin,low,high)
    topadded  = TH1D('topadded', 'topadded', bin, low, high)
    wjetsadded = TH1D('wjetsadded', 'wjetsadded',bin,low,high)
    wjetsaddedpre = TH1D('wjetsaddedpre', 'wjetsaddedpre',bin,low,high)

    wjetsadded = Variables['WJets'].Clone()
    wjetsadded.Add(Variables['WW'])
    wjetsadded.Add(Variables['ZJets_M50'])

    wjetsaddedpre = VariablesPre['WJets'].Clone()
    wjetsaddedpre.Add(VariablesPre['WW'])
    wjetsaddedpre.Add(VariablesPre['ZJets_M50'])

    topadded = Variables['TTbar_Madgraph'].Clone()
    topadded.Add(Variables['T_s'])
    topadded.Add(Variables['Tbar_s'])
    topadded.Add(Variables['T_t'])
    topadded.Add(Variables['Tbar_t'])
    topadded.Add(Variables['T_tW'])
    topadded.Add(Variables['Tbar_tW'])

    scaledwjets  = TH1D('scaledwjets', 'scaledwjets', bin, low, high)    
    scaledwjets = VariablesPre['WJets'].Clone()
    if (VariablesPre['WJets'].Integral()>0):
        scaledwjets.Scale(Variables['WJets'].Integral()/VariablesPre['WJets'].Integral())
    scaledzjets  = TH1D('scaledzjets', 'scaledzjets', bin, low, high)    
    scaledzjets = VariablesPre['ZJets_M50'].Clone()
    if (VariablesPre['ZJets_M50'].Integral()>0):
        scaledzjets.Scale(Variables['ZJets_M50'].Integral()/VariablesPre['ZJets_M50'].Integral())
    scaledww  = TH1D('scaledww', 'scaledww', bin, low, high)    
    scaledww = VariablesPre['WW'].Clone()
    if (VariablesPre['WW'].Integral()>0):
        scaledww.Scale(Variables['WW'].Integral()/VariablesPre['WW'].Integral())

    wjetsaddedpre.Scale(wjetsadded.Integral()/wjetsaddedpre.Integral())
    wjetsaddedpre.SetFillColor(ROOT.kGreen-3)
    wjetsaddedpre.SetLineColor(1)
    wjetsaddedpre.SetLineWidth(2)
    topadded.SetFillColor(ROOT.kRed-7)
    topadded.SetLineColor(1)
    topadded.SetLineWidth(2)

    #stack.Add(wjetsaddedpre)
    #stack.Add(topadded)
    #added.Add(wjetsaddedpre)  
    #added.Add(topadded)
    stack.Add(wjetsadded)
    stack.Add(topadded)
    added.Add(wjetsadded)  
    added.Add(topadded)

    if (channel == 'electron'):
        print 'Data: ',Variables['Data_el'].Integral(),' Background: ',added.Integral(),' Data/Background: ',Variables['Data_el'].Integral()/added.Integral()
    if (channel == 'muon'):
        print 'Data: ',Variables['Data_mu'].Integral(),' Background: ',added.Integral(),' Data/Background: ',Variables['Data_mu'].Integral()/added.Integral()

    lumi_error = 0.045
    ttbar_error = 0.15
    wjets_error = 0.20
    other_error = 0.20

    uncert_list = []
    lumiband = added.Clone();
    
    List_1 = ['WW', 'ZJets_M50', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']

    for hbin in range(0,lumiband.GetNbinsX()+1): 

        uncert_lumi = 0 
        uncert_xsec = 0
        uncert_stat = 0

        for i in List_1:

            error = 0
             
            if i in xsec.keys():
                if i.startswith('T'): 
                    error = ttbar_error
                    uncert_lumi += (efficiency[i]*xsec[i])**2 
                    uncert_xsec += (Variables[i].GetBinContent(hbin+1)*error)**2
                    uncert_stat += Variables[i].GetBinError(hbin+1)**2
                elif i=='WJets':
                    error = wjets_error
                    uncert_lumi += (efficiency[i]*xsec[i])**2 
                    uncert_xsec += (Variables[i].GetBinContent(hbin+1)*error)**2
                    uncert_stat += Variables[i].GetBinError(hbin+1)**2

        uncert = sqrt( (lumi_error**2)*uncert_lumi + uncert_xsec + uncert_stat )

        if lumiband.GetBinContent(hbin+1) != 0:
            dummy = 1.0
            #print lumiband.GetBinContent(hbin+1),'+/-',uncert,'(',100*uncert/lumiband.GetBinContent(hbin+1),'%)'
        lumiband.SetBinError(hbin+1,uncert);
        added.SetBinError(hbin+1,uncert);
        uncert_list . append(uncert)
    
    
    lumiband.SetFillStyle(3344);
    lumiband.SetFillColor(1);
    gStyle.SetHatchesLineWidth(1)
        
    legend = TLegend(.60,.70,.90,.90)
    if (channel == 'electron'):
        legend . AddEntry( Variables['Data_el'], 'Data' , "lp")
    if (channel == 'muon'):
        legend . AddEntry( Variables['Data_mu'], 'Data' , "lp")

    legend . AddEntry( Variables['TTbar_Madgraph'], "t#bar{t} + Single-Top", "f")
    legend . AddEntry( Variables['WJets'], "W#rightarrowl#nu + Z/#gamma*#rightarrowl^{+}l^{-} + VV" , "f")
    #legend . AddEntry( Variables['Wprime1700Right'], "W'_{R} x 20, m=1.7 TeV", "l")
    #legend . AddEntry( Variables['Wprime1900Right'], "W'_{R} x 20, m=1.9 TeV", "l")
    #legend . AddEntry( Variables['Wprime2100Right'], "W'_{R} x 20, m=2.1 TeV", "l")
    legend . AddEntry( lumiband , "Uncertainty" , "f")

    c4 = TCanvas("c4","c4", 1000, 800)
    
    c4.SetBottomMargin(0.3)
    c4.SetRightMargin(0.06)
    stack.SetMaximum( 2*stack.GetMaximum() ) 
    if setLog:
        c4.SetLogy()
        stack.SetMaximum( stack.GetMaximum()  +  10*stack.GetMaximum() ) 

    stack.SetMinimum(0.1 ) 
    stack.Draw()
    stack.GetYaxis().CenterTitle()
    stack.GetYaxis().SetTitle(ylabel)
    stack.GetXaxis().SetLabelSize(0)
    #stack.GetXaxis().SetTitle(xlabel)
    lumiband.Draw("samee2")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")    

    if (channel == 'electron'):
        if (Variables['Data_el'].GetBinContent(bin+1)>0):
            print "Overflow for electron data in ",var
        Variables['Data_el'].SetMarkerStyle(20)
        Variables['Data_el'].Draw('SAMES:E1')
    if (channel == 'muon'):
        if (Variables['Data_mu'].GetBinContent(bin+1)>0):
            print "Overflow for muon data in ",var
        Variables['Data_mu'].SetMarkerStyle(20)
        Variables['Data_mu'].Draw('SAMES:E1')

    #Variables['Wprime1700Right'].Draw("same")
    #Variables['Wprime1900Right'].Draw("same")
    #Variables['Wprime2100Right'].Draw("same")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextAlign(31) # align right
    if channel == ('electron'): 
        latex2.DrawLatex(0.87, 0.95, "CMS Preliminary, "+lumiPlot_el+" fb^{-1} at #sqrt{s} = 8 TeV");
    if channel == ('muon'): 
        latex2.DrawLatex(0.87, 0.95, "CMS Preliminary, "+lumiPlot_mu+" fb^{-1} at #sqrt{s} = 8 TeV");

    latex3 = TLatex()
    latex3.SetNDC()
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31) # align right

    if (channel == 'electron'):
        latex3.DrawLatex(0.47, 0.85, "e+jets " + njets);   
    if (channel == 'muon'):
        latex3.DrawLatex(0.47, 0.85, "#mu+jets " + njets);   
 
    latex4 = TLatex()
    latex4.SetNDC()
    latex4.SetTextSize(0.03)
    latex4.SetTextAlign(22) # align right
    latex4.SetTextAngle(90) # align right
    #latex4.DrawLatex(0.905, 0.56, "overflow")
    
    Pull  = TH1D('Pull', 'Pull', bin, low, high)  
    if (channel == 'electron'):
        Pull = Variables['Data_el'].Clone();
    if (channel == 'muon'):
        Pull = Variables['Data_mu'].Clone();    
    Pull.Add(added,-1)
    #Pull.Divide(added)

    #for i in range(bin):
    #    i += 1
    #    #print i+1,' ',added.GetBinContent(i+1)
    #    if Pull.GetBinContent(i+1) != 0:
    #        Pull.SetBinContent(i+1,Pull.GetBinContent(i+1)/Pull.GetBinError(i+1))
    #    else: Pull.SetBinContent(i+1,0)
    for i in range(bin):
        i += 1
        #print i+1,' ',added.GetBinContent(i+1)
        if Pull.GetBinContent(i) != 0:
            Pull.SetBinContent(i,Pull.GetBinContent(i)/Pull.GetBinError(i))
        else: Pull.SetBinContent(i,0)

    pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
    pad.SetTopMargin(0.7)
    pad.SetFillColor(0)
    pad.SetGridy(1)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)
    pad.SetRightMargin(0.06)

    Pull.SetMarkerStyle(20)
    Pull.SetMaximum(3.0 )
    Pull.SetMinimum(-3.0)
    Pull.SetFillColor(2)
    Pull.GetXaxis().SetTitle(xlabel)
    Pull.GetYaxis().SetTitleSize(0.04)
    Pull.GetYaxis().SetTitle('#sigma(Data-MC)')
    Pull.SetMarkerSize(0.7)
    Pull.GetYaxis().SetNdivisions(5);
    Pull.Draw("HIST")
    
    if (finalcuts): save = save+"_finalcuts"

    c4.SaveAs('Wprime_Plots/StackedHisto_' + save + '.pdf')
    #c4.SaveAs('Wprime_Plots/StackedHisto_' + save + '.png')
    #c4.SaveAs('Wprime_Plots/StackedHisto_' + save + '.C')
      
    del c4
    #del pad

    
    rebn = 2
    added.Rebin(rebn)

    if (channel == 'electron'):
        Variables['Data_el'].Rebin(rebn)
        KSN = Variables["Data_el"].KolmogorovTest(added,"ON")
        KS = Variables["Data_el"].KolmogorovTest(added,"OD")
        KSM = Variables["Data_el"].KolmogorovTest(added,"OM")
        ch2 = Variables['Data_el'].Chi2Test(added, "UWP")

    if (channel == 'muon'):
        Variables['Data_mu'].Rebin(rebn)
        KSN = Variables["Data_mu"].KolmogorovTest(added,"ON")
        KS = Variables["Data_mu"].KolmogorovTest(added,"OD")
        KSM = Variables["Data_mu"].KolmogorovTest(added,"OM")
        ch2 = Variables['Data_mu'].Chi2Test(added, "UWP")

    ksr = round(KS, 4)
    ksrn = round(KSN, 4)
    ksrm = round(KSM, 4)
    ch2r = round(ch2, 4)

    print 'var =', var, ' KS =', ksr
    print 'var =', var, ' Chi2 =', ch2r

    strks = save + ' KS = ' + str(ksr) + ' KSN = ' + str(ksrn) + ' KSM = ' + str(ksrm) + ' CHI2pvalue = ' + str(ch2r)
    ksfile.write(strks+"\n")

    del stack
    del added
    del Variables



######################
## electron channel ##
######################
channel = 'electron'

var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_el'; setLog = True; nBtags = 2; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 2; finalcuts = True;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 3;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'Jet1Jet2W_M_el'; setLog = True; nBtags = 2; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 3;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'corr_met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'corr_met_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'type1corrmet_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'type1corrmet_el'; setLog = True; nBtags = 2; finalcuts = False;    
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'met_el'; setLog = True; nBtags = 2; finalcuts = False;    
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_0_pt_WprimeCalc'; bin = 130; low = 0; high = 1300; xaxis = "p_{T} (jet1) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet0_pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_0_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet1)"; yaxis = 'Events / 0.1'; save = 'PFjet0_eta_el'; setLog = False; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_1_pt_WprimeCalc'; bin = 120; low = 0; high = 1200; xaxis = "p_{T} (jet2) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet1_pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_1_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet2)"; yaxis = 'Events / 0.1'; save = 'PFjet1_eta_el'; setLog = False; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_2_pt_WprimeCalc'; bin = 80; low = 0; high = 800; xaxis = "p_{T} (jet3) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet2_pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_2_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet3)"; yaxis = 'Events / 0.1'; save = 'PFjet2_eta_el'; setLog = False; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'elec_1_pt_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "electron p_{T} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'electron0_pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'elec_1_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "electron #eta"; yaxis = 'Events / 0.1'; save = 'electron0_eta_el'; setLog = False; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'elec_1_RelIso_WprimeCalc'; bin =30; low = 0; high = 0.15; xaxis = "PF Rel. Isolation"; yaxis = 'Events / 0.005'; save = 'relIso_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_LjetsTopoCalcNew'; bin = 100; low = 0; high = 1000; xaxis = "M(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopMass_Best_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_Pt_LjetsTopoCalcNew'; bin = 150; low = 0; high = 1500; xaxis = "p_{T}(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopPt_Best_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2_Pt_LjetsTopoCalcNew'; bin = 150; low = 0; high = 1500; xaxis = "p_{T}(jet1,jet2) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'Pt_Jet1Jet2_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Ht_LjetsTopoCalcNew'; bin = 125; low = 0; high = 2500; xaxis = "HT [GeV]"; yaxis = 'Events / 20 GeV'; save = 'HT_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'nPV_WprimeCalc'; bin = 50; low = 0; high = 50; xaxis = "# Vertices"; yaxis = 'Events'; save = 'nPV_el'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'n_btags_WprimeCalc'; bin = 8; low = -0.5; high = 7.5; xaxis = "# B-tags"; yaxis = 'Events'; save = 'nbtags_el'; setLog = True; nBtags = -1; finalcuts = False;    
plot_forPAS(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'dphiLepJ1_LjetsTopoCalcNew'; bin = 64; low = 0.0; high = 3.2; xaxis = "#Delta#phi(el,jet1)"; yaxis = 'Events'; save = 'dphiLepJ1_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'DphiJMET_LjetsTopoCalcNew'; bin = 64; low = 0.0; high = 3.2; xaxis = "#Delta#phi(MET,jet1)"; yaxis = 'Events'; save = 'DphJMET_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_3_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet4)"; yaxis = 'Events'; save = 'PFjet3_eta_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_2_pt_WprimeCalc'; bin = 60; low = 0; high = 600; xaxis = "p_{T} (jet3) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet2_pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'AllJets_M_LjetsTopoCalcNew'; bin = 100; low = 0; high = 5000; xaxis = "M(all jets) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'AllJets_M_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'AplanarityMu_LjetsTopoCalcNew'; bin = 100; low = 0; high = 0.5; xaxis = "Aplanarity(all jets,el)"; yaxis = 'Events / 0.005'; save = 'AplanarityMu_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Sphericity_LjetsTopoCalcNew'; bin = 100; low = 0; high = 1.0; xaxis = "Sphericity(all jets)"; yaxis = 'Events / 0.01'; save = 'Sphericity_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1TagJet2TagW_M_LjetsTopoCalcNew'; bin = 90; low = 0; high = 4500; xaxis = "M(jet1 tag, jet2 tag,W) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'Jet1TagJet2TagW_M_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2_M_LjetsTopoCalcNew'; bin = 80; low = 0; high = 4000; xaxis = "M(jet1, jet2) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'Jet1Jet2_M_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'SecBTagTop_Pt_LjetsTopoCalcNew'; bin = 62; low = -50.0; high = 1500; xaxis = "Pt(2nd btag top) [GeV]"; yaxis = 'Events / 25 GeV'; save = 'SecBTagTop_Pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'SecBTagTopMass_LjetsTopoCalcNew'; bin = 51; low = -50.0; high = 2500; xaxis = "Pt(2nd btag top) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'SecBTagTop_Pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

#BAD
var = 'LeptonJet_DeltaR_LjetsTopoCalcNew'; bin = 60; low = 0.0; high = 6.0; xaxis = "#DeltaR(el,jet1 or jet2)"; yaxis = 'Events / 0.1'; save = 'LeptonJet_DeltaR_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

#BAD
var = 'J2_NotBestJet_Pt_LjetsTopoCalcNew'; bin = 30; low = -200.0; high = 1300; xaxis = "Pt(jet2 not best) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'J2_NotBestJet_Pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BTagTop_Pt_LjetsTopoCalcNew'; bin = 62; low = -50.0; high = 1500; xaxis = "Pt(btag top) [GeV]"; yaxis = 'Events / 25 GeV'; save = 'BTagTop_Pt_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'SqrtsT_LjetsTopoCalcNew'; bin = 60; low = 0.0; high = 3000; xaxis = "sqrt(sT) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'SqrtsT_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'PzOverHT_LjetsTopoCalcNew'; bin = 60; low = -6.0; high = 6.0; xaxis = "pz/HT"; yaxis = 'Events / 50 GeV'; save = 'PzOverHT_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2_DeltaPhi_LjetsTopoCalcNew'; bin = 64; low = 0.0; high = 3.2; xaxis = "#Delta#phi(jet1,jet2)"; yaxis = 'Events / 0.05'; save = 'Jet1Jet2_DeltaPhi_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Cos_BestJetLepton_BestTop_LjetsTopoCalcNew'; bin = 100; low = -1.0; high = 1.0; xaxis = "cos(bestjet,el)"; yaxis = 'Events / 0.02'; save = 'Cos_BestJetLepton_BestTop_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

#BAD
var = 'Cos_LightjetJetLepton_BestTop_LjetsTopoCalcNew'; bin = 650; low = -12.0; high = 1.0; xaxis = "cos(lightjet,el)"; yaxis = 'Events / 0.02'; save = 'Cos_BestJetLepton_BestTop_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)


######################
## muon channel ##
######################
channel = 'muon'

var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_mu'; setLog = True; nBtags = 2; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 2; finalcuts = True;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 3;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'Jet1Jet2W_M_mu'; setLog = True; nBtags = 2; finalcuts = False;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = 3;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'corr_met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'corr_met_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'type1corrmet_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'type1corrmet_mu'; setLog = True; nBtags = 2; finalcuts = False;    
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'met_mu'; setLog = True; nBtags = 2; finalcuts = False;    
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
##plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_0_pt_WprimeCalc'; bin = 130; low = 0; high = 1300; xaxis = "p_{T} (jet1) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet0_pt_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_0_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet1)"; yaxis = 'Events / 0.1'; save = 'PFjet0_eta_mu'; setLog = False; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_1_pt_WprimeCalc'; bin = 120; low = 0; high = 1200; xaxis = "p_{T} (jet2) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet1_pt_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_1_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet2)"; yaxis = 'Events / 0.1'; save = 'PFjet1_eta_mu'; setLog = False; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_2_pt_WprimeCalc'; bin = 80; low = 0; high = 800; xaxis = "p_{T} (jet3) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'PFjet2_pt_mu'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'jet_2_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "#eta (jet3)"; yaxis = 'Events / 0.1'; save = 'PFjet2_eta_mu'; setLog = False; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'muon_1_pt_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "muon p_{T} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'muon_pt_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'muon_1_eta_WprimeCalc'; bin = 50; low = -2.5; high = 2.5; xaxis = "muon #eta"; yaxis = 'Events / 0.1'; save = 'muon_eta_mu'; setLog = False; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'muon_1_RelIso_WprimeCalc'; bin =30; low = 0; high = 0.15; xaxis = "PF Rel. Isolation"; yaxis = 'Events / 0.005'; save = 'relIso_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_LjetsTopoCalcNew'; bin = 100; low = 0; high = 1000; xaxis = "M(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopMass_Best_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_Pt_LjetsTopoCalcNew'; bin = 150; low = 0; high = 1500; xaxis = "p_{T}(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopPt_Best_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Jet1Jet2_Pt_LjetsTopoCalcNew'; bin = 150; low = 0; high = 1500; xaxis = "p_{T}(jet1,jet2) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'Pt_Jet1Jet2_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
NBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'Ht_LjetsTopoCalcNew'; bin = 125; low = 0; high = 2500; xaxis = "HT [GeV]"; yaxis = 'Events / 20 GeV'; save = 'HT_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'nPV_WprimeCalc'; bin = 50; low = 0; high = 50; xaxis = "# Vertices"; yaxis = 'Events'; save = 'nPV_mu'; setLog = True; nBtags = 2; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'n_btags_WprimeCalc'; bin = 8; low = -0.5; high = 7.5; xaxis = "# B-tags"; yaxis = 'Events'; save = 'nbtags_mu'; setLog = True; nBtags = -1; finalcuts = False;    
plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'LeptonJet_DeltaR_LjetsTopoCalcNew'; bin = 60; low = 0.0; high = 6.0; xaxis = "#DeltaR(el,jet1 or jet2)"; yaxis = 'Events / 0.1'; save = 'LeptonJet_DeltaR_mu'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_forPAS(channel,var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

ksfile.close()


def plot_CompareTTbar(channel, var, bin, low, high, ylabel, xlabel, save, nBtags = 2, setLog = False, finalcuts = False):

    List = ['TTbar_Madgraph','TTbar_Powheg']
       
    Variables = {}
    
    njets = ""

    if nBtags == -1:
        cutbtag = ''
        save = save     
    if nBtags == 0:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) == 0 )'
        save = save + '_0bTags'
        njets = "N_{b tags} = 0"
    if nBtags == 1:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) == 1)'
        save = save + '_1bTags'
        njets = "N_{b tags} = 1"
    if nBtags == 2:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) >= 1 ) '
        save = save + '_GE1bTags'
        njets = "N_{b tags} #geq 1"


    legend = TLegend(.70,.70,.90,.90)

    background = 0
    efficiency = {}

    j = 0
    k = 0
    for Type in List:
        Variables[Type] = TH1D(Type+var+channel, Type+var+channel, bin, low, high)

        if (channel == 'electron'):      
            cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 30 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 
        if (channel == 'muon'):      
            cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20' 

        if finalcuts: cut+= '&& BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 75  && Jet1Jet2_Pt_LjetsTopoCalcNew > 100'  

        histName = Type+var+channel

        Trees[Type].Draw(var + " >> " + histName, "(weight_PU_PileUpCalc)*(" + cut + cutbtag + ")", 'goff')
            
        if Variables[Type].Integral() != 0:
            Variables[Type].Scale ( 1/Variables[Type].Integral() ) 

    Variables['TTbar_Madgraph'].SetLineColor(2)
    Variables['TTbar_Madgraph'].SetLineWidth(3)

    Variables['TTbar_Powheg'].SetLineColor(4)
    Variables['TTbar_Powheg'].SetLineWidth(3)


    legend . AddEntry( Variables['TTbar_Madgraph'], 't#bar{t} Madgraph', "l")
    legend . AddEntry( Variables['TTbar_Powheg'], 't#bar{t} Powheg' , "l")

 
    c4 = TCanvas("c4","Full PF Cuts for Data vs WJets MC", 1000, 800)
    
    c4.SetRightMargin(0.06)
    if setLog:
        c4.SetLogy()

    Variables['TTbar_Madgraph'].GetXaxis().SetTitle(xlabel)
    Variables['TTbar_Madgraph'].GetYaxis().SetTitle("arbitrary units")

    Variables['TTbar_Madgraph'].SetMaximum(Variables['TTbar_Madgraph'].GetMaximum() + 10*Variables['TTbar_Madgraph'].GetMaximum())
    Variables['TTbar_Madgraph'].Draw()
    Variables['TTbar_Powheg'].Draw("same")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.87, 0.95, "CMS Simulation");

    latex3 = TLatex()
    latex3.SetNDC()
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31) # align right
    if (channel == 'electron'):
        latex3.DrawLatex(0.55, 0.85, "e+jets " + njets);   
    if (channel == 'muon'):
        latex3.DrawLatex(0.55, 0.85, "#mu+jets " + njets);   
 
    c4.SaveAs('Wprime_Plots/compareTTbar_' + save + '.pdf')
    c4.SaveAs('Wprime_Plots/compareTTbar_' + save + '.png')
     
    del c4
    #del pad
    del Variables

channel = 'electron'
var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_el'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_LjetsTopoCalcNew'; bin = 100; low = 0; high = 1000; xaxis = "M(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopMass_Best_el'; setLog = True; nBtags = 2; finalcuts = False;   
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'n_btags_WprimeCalc'; bin = 8; low = -0.5; high = 7.5; xaxis = "# B-tags"; yaxis = 'Events'; save = 'nbtags_el'; setLog = True; nBtags = -1; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'nSelJets_CommonCalc'; bin = 15; low = -0.5; high = 14.5; xaxis = "# Jets"; yaxis = 'Events'; save = 'njets_el'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

channel = 'muon'
var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_mu'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'BestTop_LjetsTopoCalcNew'; bin = 100; low = 0; high = 1000; xaxis = "M(best jet,W) [GeV]"; yaxis = 'Events / 10 GeV'; save = 'TopMass_Best_mu'; setLog = True; nBtags = 2; finalcuts = False;   
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'n_btags_WprimeCalc'; bin = 8; low = -0.5; high = 7.5; xaxis = "# B-tags"; yaxis = 'Events'; save = 'nbtags_mu'; setLog = True; nBtags = -1; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'nSelJets_CommonCalc'; bin = 15; low = -0.5; high = 14.5; xaxis = "# Jets"; yaxis = 'Events'; save = 'njets_mu'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareTTbar(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

def plot_CompareWJetsHF(channel, var, bin, low, high, ylabel, xlabel, save, nBtags = 2, setLog = False, finalcuts = False):
       
    Variables = {}
    
    njets = ""

    if nBtags == -1:
        cutbtag = ''
        save = save     
    if nBtags == 0:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) == 0 )'
        save = save + '_0bTags'
        njets = "N_{b tags} = 0"
    if nBtags == 1:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) == 1)'
        save = save + '_1bTags'
        njets = "N_{b tags} = 1"
    if nBtags == 2:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) >= 1 ) '
        save = save + '_GE1bTags'
        njets = "N_{b tags} #geq 1"

    cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
    cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
    cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light
 
    legend = TLegend(.70,.70,.90,.90)

    background = 0
    efficiency = {}

    WjjHist = TH1D('WjjHist'+var+channel, 'WjjHist'+var+channel, bin, low, high)
    WccHist = TH1D('WccHist'+var+channel, 'WccHist'+var+channel, bin, low, high)
    WbbHist = TH1D('WbbHist'+var+channel, 'WbbHist'+var+channel, bin, low, high)

    if (channel == 'electron'):      
        cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 32 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 
    if (channel == 'muon'):      
        cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20' 

    if finalcuts: cut+= '&& BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 75  && Jet1Jet2_Pt_LjetsTopoCalcNew > 100'  

    if (channel=='electron'):
        weight = 'weight_PU_ABC_PileUpCalc*weight_ElectronEff_WprimeCalc'
    if (channel=='muon'):
        weight = 'weight_PU_ABC_PileUpCalc*weight_MuonEff_WprimeCalc'


    Trees['WJets'].Draw(var + " >> " + 'WjjHist'+var+channel, "("+weight+")*(" + cut + cutbtag + cutwjj + ")", 'goff')
    Trees['WJets'].Draw(var + " >> " + 'WbbHist'+var+channel, "("+weight+")*(" + cut + cutbtag + cutwbb + ")", 'goff')
    Trees['WJets'].Draw(var + " >> " + 'WccHist'+var+channel, "("+weight+")*(" + cut + cutbtag + cutwcc + ")", 'goff')            

    if WjjHist.Integral() != 0:
        WjjHist.Scale ( 1/WjjHist.Integral() ) 
    if WbbHist.Integral() != 0:
        WbbHist.Scale ( 1/WbbHist.Integral() ) 
    if WccHist.Integral() != 0:
        WccHist.Scale ( 1/WccHist.Integral() ) 


    WjjHist.SetLineColor(1)
    WjjHist.SetLineWidth(2)

    WbbHist.SetLineColor(2)
    WbbHist.SetLineWidth(2)

    WccHist.SetLineColor(4)
    WccHist.SetLineWidth(2)


    legend . AddEntry( WjjHist, 'W + light-jets', "l")
    legend . AddEntry( WccHist, 'W + c-jets' , "l")
    legend . AddEntry( WbbHist, 'W + b-jets' , "l")

 
    c4 = TCanvas("c4","Full PF Cuts for Data vs WJets MC", 1000, 800)
    
    c4.SetRightMargin(0.06)
    if setLog:
        c4.SetLogy()

    WjjHist.GetXaxis().SetTitle(xlabel)
    WjjHist.GetYaxis().SetTitle("arbitrary units")

    WjjHist.SetMaximum(WjjHist.GetMaximum() + 10*WjjHist.GetMaximum())
    WjjHist.Draw()
    WccHist.Draw("same")
    WbbHist.Draw("same")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.87, 0.95, "CMS Simulation");

    latex3 = TLatex()
    latex3.SetNDC()
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31) # align right
    if (channel == 'electron'):
        latex3.DrawLatex(0.55, 0.85, "e+jets " + njets);   
    if (channel == 'muon'):
        latex3.DrawLatex(0.55, 0.85, "#mu+jets " + njets);   
 
    c4.SaveAs('Wprime_Plots/compareWJetsHF_' + save + '.pdf')
    c4.SaveAs('Wprime_Plots/compareWJetsHF_' + save + '.png')
     
    del c4
    #del pad
    del Variables

###################
### electrons #####
###################
channel = 'electron'

var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_el'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'elec_1_RelIso_WprimeCalc'; bin =30; low = 0; high = 0.15; xaxis = "PF Rel. Isolation"; yaxis = 'Events / 0.005'; save = 'relIso_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'corr_met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'corr_met_el'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)


###################
##### muons #######
###################
channel = 'muon'

var = 'BestJetJet2W_M_LjetsTopoCalcNew'; bin = 68; low = 100; high = 3500; xaxis = "M(tb) [GeV]"; yaxis = 'Events / 50 GeV'; save = 'BestJetJet2W_M_mu'; setLog = True; nBtags = 2; finalcuts = False;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'muon_1_RelIso_WprimeCalc'; bin =30; low = 0; high = 0.15; xaxis = "PF Rel. Isolation"; yaxis = 'Events / 0.005'; save = 'relIso_mu'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)

var = 'corr_met_WprimeCalc'; bin = 100; low = 0; high = 1000; xaxis = "E_{T}^{miss} [GeV]"; yaxis = 'Events / 10 GeV'; save = 'corr_met_mu'; setLog = True; nBtags = 2; finalcuts = False;    
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)
nBtags = -1;
#plot_CompareWJetsHF(channel, var, bin, low, high, yaxis, xaxis , save, nBtags, setLog, finalcuts)



def plot_data2D(channel, var1, var2, xbin, xlow, xhigh, ybin, ylow, yhigh, xlabel, ylabel, nBtags = 2, finalcuts = False):

    if (channel == 'electron'): List = ['Data_el']
    if (channel == 'muon'): List = ['Data_mu']
       
    Variables = {}
    
    njets = ""

    if nBtags == -1:
        cutbtag = ''
        save = ''     
    if nBtags == 0:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) == 0 )'
        save = '_0bTags'
        njets = "N_{b tags} = 0"
    if nBtags == 1:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) == 1)'
        save = '_1bTags'
        njets = "N_{b tags} = 1"
    if nBtags == 2:
        cutbtag =  ' && ( (jet_0_tag_WprimeCalc + jet_1_tag_WprimeCalc) >= 1 ) '
        save = '_GE1bTags'
        njets = "N_{b tags} #geq 1"



    background = 0
    efficiency = {}

    j = 0
    k = 0
    for Type in List:
        Variables[Type] = TH2D(Type+var+channel, Type+var+channel,xbin,xlow,xhigh,ybin,ylow,yhigh)

        if (channel == 'electron'):      
            cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 30 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 
        if (channel == 'muon'):      
            cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20' 

        if finalcuts: cut+= '&& BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 75  && Jet1Jet2_Pt_LjetsTopoCalcNew > 100'  

        histName = Type+var+channel
        Trees[Type].Draw(var2+":"+var1 + " >> " + histName, "(" + cut + cutbtag + ")", 'goff')
            

    if (channel == 'electron'): 
        chan = '_el'
        Datatype = 'Data_el'
    if (channel == 'muon'): 
        Datatype = 'Data_mu'
        chan = '_mu'

    print Variables[Datatype].Integral()
    Variables[Datatype].SetFillColor(1)

    legend = TLegend(.70,.80,.90,.90)
    legend . AddEntry( Variables[Datatype], 'Data', "l")

 
    c4 = TCanvas("c4","Full PF Cuts for Data vs WJets MC", 1000, 800)
    c4.SetRightMargin(0.06)

    Variables[Datatype].GetXaxis().SetTitle(xlabel)
    Variables[Datatype].GetYaxis().SetTitle(ylabel)

    #Variables[Datatype].Draw("box")
    Variables[Datatype].Draw("scat")
    #Variables[Datatype].Draw("colz")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.87, 0.95, "CMS Preliminary, "+lumiPlot_el+" fb^{-1} at #sqrt{s} = 8 TeV");

    latex3 = TLatex()
    latex3.SetNDC()
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31) # align right
    if (channel == 'electron'):
        latex3.DrawLatex(0.55, 0.85, "e+jets " + njets);   
    if (channel == 'muon'):
        latex3.DrawLatex(0.55, 0.85, "#mu+jets " + njets);   
 
    c4.SaveAs('Wprime_Plots/data2D_' + var1 + 'vs' + var2 + save + chan + '.pdf')
    c4.SaveAs('Wprime_Plots/data2D_' + var1 + 'vs' + var2 + save + chan + '.png')
     
    del c4
    #del pad
    del Variables

channel = 'electron'; var1 = 'corr_met_WprimeCalc'; var2 = 'elec_1_RelIso_WprimeCalc';
xbin = 100; xlow = 0; xhigh = 1000; xlabel = "E_{T}^{miss} [GeV]";
ybin = 30; ylow = 0; yhigh = 0.15; ylabel = "PF Rel. Isolation";
#plot_data2D(channel, var1, var2, xbin, xlow, xhigh, ybin, ylow, yhigh, xlabel, ylabel, nBtags = 2, finalcuts = False)

channel = 'muon'; var1 = 'corr_met_WprimeCalc'; var2 = 'muon_1_RelIso_WprimeCalc';
xbin = 100; xlow = 0; xhigh = 1000; xlabel = "E_{T}^{miss} [GeV]";
ybin = 30; ylow = 0; yhigh = 0.15; ylabel = "PF Rel. Isolation";
#plot_data2D(channel, var1, var2, xbin, xlow, xhigh, ybin, ylow, yhigh, xlabel, ylabel, nBtags = 2, finalcuts = False)





global syncfile
syncfile = open('myevents.txt', 'w')
global highmassfile
highmassfile = open('highmassevents.txt','w')

def Synch():
    
    cnt = 0
    zerobtags = 0
    onebtags = 0
    ge1btags = 0
    ge2btags = 0
    final = 0

    List = ['TTbar_Madgraph']
    for Type in List:
    
        if Type == 'TTbar_Madgraph':

            Nevent = Trees[Type].GetEntriesFast()        
          
            for entry in xrange(Nevent):
                Trees[Type].GetEntry(entry)            

                event = Trees[Type].event_CommonCalc
                run = Trees[Type].run_CommonCalc
                lumisec = Trees[Type].lumi_CommonCalc

                jet0_pt=Trees[Type].jet_0_pt_WprimeCalc
                jet0_eta=Trees[Type].jet_0_eta_WprimeCalc
                jet0_phi=Trees[Type].jet_0_phi_WprimeCalc
                jet1_pt=Trees[Type].jet_1_pt_WprimeCalc
                jet1_eta=Trees[Type].jet_1_eta_WprimeCalc
                jet1_phi=Trees[Type].jet_1_phi_WprimeCalc
                jet2_pt=Trees[Type].jet_2_pt_WprimeCalc
                jet2_eta=Trees[Type].jet_2_eta_WprimeCalc
                jet2_phi=Trees[Type].jet_2_phi_WprimeCalc
                jet3_pt=Trees[Type].jet_3_pt_WprimeCalc
                jet3_eta=Trees[Type].jet_3_eta_WprimeCalc
                jet3_phi=Trees[Type].jet_3_phi_WprimeCalc
                jet4_pt=Trees[Type].jet_4_pt_WprimeCalc
                jet4_eta=Trees[Type].jet_4_eta_WprimeCalc
                jet4_phi=Trees[Type].jet_4_phi_WprimeCalc
                jet5_pt=Trees[Type].jet_5_pt_WprimeCalc
                jet5_eta=Trees[Type].jet_5_eta_WprimeCalc
                jet5_phi=Trees[Type].jet_5_phi_WprimeCalc
                jet6_pt=Trees[Type].jet_6_pt_WprimeCalc
                jet6_eta=Trees[Type].jet_6_eta_WprimeCalc
                jet6_phi=Trees[Type].jet_6_phi_WprimeCalc
                jet7_pt=Trees[Type].jet_7_pt_WprimeCalc
                jet7_eta=Trees[Type].jet_7_eta_WprimeCalc
                jet7_phi=Trees[Type].jet_7_phi_WprimeCalc
                jet8_pt=Trees[Type].jet_8_pt_WprimeCalc
                jet8_eta=Trees[Type].jet_8_eta_WprimeCalc
                jet8_phi=Trees[Type].jet_8_phi_WprimeCalc
                jet9_pt=Trees[Type].jet_9_pt_WprimeCalc
                jet9_eta=Trees[Type].jet_9_eta_WprimeCalc
                jet9_phi=Trees[Type].jet_9_phi_WprimeCalc
 
                btag0=Trees[Type].jet_0_tag_WprimeCalc
                btag1=Trees[Type].jet_1_tag_WprimeCalc

                muon_pt = Trees[Type].muon_1_pt_WprimeCalc
                muon_eta = Trees[Type].muon_1_eta_WprimeCalc
                muon_phi = Trees[Type].muon_1_phi_WprimeCalc

                electron_pt = Trees[Type].elec_1_pt_WprimeCalc
                electron_eta = Trees[Type].elec_1_eta_WprimeCalc
                electron_phi = Trees[Type].elec_1_phi_WprimeCalc

                corr_met = Trees[Type].corr_met_WprimeCalc
                met = Trees[Type].met_WprimeCalc
 
                Jet1Jet2W_M=Trees[Type].Jet1Jet2W_M_LjetsTopoCalcNew
                BestTop_M=Trees[Type].BestTop_LjetsTopoCalcNew

                 
                muon_relIso = Trees[Type].muon_1_RelIso_WprimeCalc
                electron_relIso = Trees[Type].elec_1_RelIso_WprimeCalc

                nbtags = 0


                deta0 = electron_eta - jet0_eta
                dphi0 = electron_phi - jet0_phi;
                if ( dphi0 >  ROOT.TMath.Pi() ): dphi0 -= 2.*ROOT.TMath.Pi();
                if ( dphi0 <= -ROOT.TMath.Pi() ): dphi0 += 2.*ROOT.TMath.Pi();
                deltaR0_el = ROOT.TMath.Sqrt(deta0*deta0 + dphi0*dphi0);      
 
                deta1 = electron_eta - jet1_eta
                dphi1 = electron_phi - jet1_phi;
                if ( dphi1 >  ROOT.TMath.Pi() ): dphi1 -= 2.*ROOT.TMath.Pi();
                if ( dphi1 <= -ROOT.TMath.Pi() ): dphi1 += 2.*ROOT.TMath.Pi();
                deltaR1_el = ROOT.TMath.Sqrt(deta1*deta1 + dphi1*dphi1);    
 
                deta2 = electron_eta - jet2_eta
                dphi2 = electron_phi - jet2_phi;
                if ( dphi2 >  ROOT.TMath.Pi() ): dphi2 -= 2.*ROOT.TMath.Pi();
                if ( dphi2 <= -ROOT.TMath.Pi() ): dphi2 += 2.*ROOT.TMath.Pi();
                deltaR2_el = ROOT.TMath.Sqrt(deta2*deta2 + dphi2*dphi2);    

                deta3 = electron_eta - jet3_eta
                dphi3 = electron_phi - jet3_phi;
                if ( dphi3 >  ROOT.TMath.Pi() ): dphi3 -= 2.*ROOT.TMath.Pi();
                if ( dphi3 <= -ROOT.TMath.Pi() ): dphi3 += 2.*ROOT.TMath.Pi();
                deltaR3_el = ROOT.TMath.Sqrt(deta3*deta3 + dphi3*dphi3);    

                deta4 = electron_eta - jet4_eta
                dphi4 = electron_phi - jet4_phi;
                if ( dphi4 >  ROOT.TMath.Pi() ): dphi4 -= 2.*ROOT.TMath.Pi();
                if ( dphi4 <= -ROOT.TMath.Pi() ): dphi4 += 2.*ROOT.TMath.Pi();
                deltaR4_el = ROOT.TMath.Sqrt(deta4*deta4 + dphi4*dphi4);    

                deta5 = electron_eta - jet5_eta
                dphi5 = electron_phi - jet5_phi;
                if ( dphi5 >  ROOT.TMath.Pi() ): dphi5 -= 2.*ROOT.TMath.Pi();
                if ( dphi5 <= -ROOT.TMath.Pi() ): dphi5 += 2.*ROOT.TMath.Pi();
                deltaR5_el = ROOT.TMath.Sqrt(deta5*deta5 + dphi5*dphi5);    

                deta6 = electron_eta - jet6_eta
                dphi6 = electron_phi - jet6_phi;
                if ( dphi6 >  ROOT.TMath.Pi() ): dphi6 -= 2.*ROOT.TMath.Pi();
                if ( dphi6 <= -ROOT.TMath.Pi() ): dphi6 += 2.*ROOT.TMath.Pi();
                deltaR6_el = ROOT.TMath.Sqrt(deta6*deta6 + dphi6*dphi6);    

                deta7 = electron_eta - jet7_eta
                dphi7 = electron_phi - jet7_phi;
                if ( dphi7 >  ROOT.TMath.Pi() ): dphi7 -= 2.*ROOT.TMath.Pi();
                if ( dphi7 <= -ROOT.TMath.Pi() ): dphi7 += 2.*ROOT.TMath.Pi();
                deltaR7_el = ROOT.TMath.Sqrt(deta7*deta7 + dphi7*dphi7);


                deta0 = muon_eta - jet0_eta
                dphi0 = muon_phi - jet0_phi;
                if ( dphi0 >  ROOT.TMath.Pi() ): dphi0 -= 2.*ROOT.TMath.Pi();
                if ( dphi0 <= -ROOT.TMath.Pi() ): dphi0 += 2.*ROOT.TMath.Pi();
                deltaR0_mu = ROOT.TMath.Sqrt(deta0*deta0 + dphi0*dphi0);    
 
                deta1 = muon_eta - jet1_eta
                dphi1 = muon_phi - jet1_phi;
                if ( dphi1 >  ROOT.TMath.Pi() ): dphi1 -= 2.*ROOT.TMath.Pi();
                if ( dphi1 <= -ROOT.TMath.Pi() ): dphi1 += 2.*ROOT.TMath.Pi();
                deltaR1_mu = ROOT.TMath.Sqrt(deta1*deta1 + dphi1*dphi1);    
 
                deta2 = muon_eta - jet2_eta
                dphi2 = muon_phi - jet2_phi;
                if ( dphi2 >  ROOT.TMath.Pi() ): dphi2 -= 2.*ROOT.TMath.Pi();
                if ( dphi2 <= -ROOT.TMath.Pi() ): dphi2 += 2.*ROOT.TMath.Pi();
                deltaR2_mu = ROOT.TMath.Sqrt(deta2*deta2 + dphi2*dphi2);    

                deta3 = muon_eta - jet3_eta
                dphi3 = muon_phi - jet3_phi;
                if ( dphi3 >  ROOT.TMath.Pi() ): dphi3 -= 2.*ROOT.TMath.Pi();
                if ( dphi3 <= -ROOT.TMath.Pi() ): dphi3 += 2.*ROOT.TMath.Pi();
                deltaR3_mu = ROOT.TMath.Sqrt(deta3*deta3 + dphi3*dphi3);    

                deta4 = muon_eta - jet4_eta
                dphi4 = muon_phi - jet4_phi;
                if ( dphi4 >  ROOT.TMath.Pi() ): dphi4 -= 2.*ROOT.TMath.Pi();
                if ( dphi4 <= -ROOT.TMath.Pi() ): dphi4 += 2.*ROOT.TMath.Pi();
                deltaR4_mu = ROOT.TMath.Sqrt(deta4*deta4 + dphi4*dphi4);    

                deta5 = muon_eta - jet5_eta
                dphi5 = muon_phi - jet5_phi;
                if ( dphi5 >  ROOT.TMath.Pi() ): dphi5 -= 2.*ROOT.TMath.Pi();
                if ( dphi5 <= -ROOT.TMath.Pi() ): dphi5 += 2.*ROOT.TMath.Pi();
                deltaR5_mu = ROOT.TMath.Sqrt(deta5*deta5 + dphi5*dphi5);    

                deta6 = muon_eta - jet6_eta
                dphi6 = muon_phi - jet6_phi;
                if ( dphi6 >  ROOT.TMath.Pi() ): dphi6 -= 2.*ROOT.TMath.Pi();
                if ( dphi6 <= -ROOT.TMath.Pi() ): dphi6 += 2.*ROOT.TMath.Pi();
                deltaR6_mu = ROOT.TMath.Sqrt(deta6*deta6 + dphi6*dphi6);    

                deta7 = muon_eta - jet7_eta
                dphi7 = muon_phi - jet7_phi;
                if ( dphi7 >  ROOT.TMath.Pi() ): dphi7 -= 2.*ROOT.TMath.Pi();
                if ( dphi7 <= -ROOT.TMath.Pi() ): dphi7 += 2.*ROOT.TMath.Pi();
                deltaR7_mu = ROOT.TMath.Sqrt(deta7*deta7 + dphi7*dphi7);

                #if (deltaR0_el < 0.3): continue
                #if (deltaR1_el < 0.3): continue
                #if (deltaR2_el < 0.3): continue
                #if (deltaR3_el < 0.3): continue
                #if (deltaR4_el < 0.3): continue
                #if (deltaR5_el < 0.3): continue
                #if (deltaR6_el < 0.3): continue
                #if (deltaR7_el < 0.3): continue
                           
                             
                #if (fabs(muon_eta) >= 2.1):
                if (fabs(electron_eta) >= 2.5):
                    #print "failed muon eta"
                    continue
                #if (muon_pt < 26):
                if (electron_pt < 30):
                    #print "failed muon pt"
                    continue
                if (jet0_pt < 100):
                    #print "failed cut 1"
                    continue
                if (jet1_pt < 40 ):
                    #print "failed cut 2"
                    continue
                #if (jet2_pt < 30 ):
                #    #print "failed cut 3"
                #    continue
                #if (jet3_pt < 30 ):
                #    #print "failed cut 4"
                #    continue
                #if (muon_relIso >= 0.12):
                if (electron_relIso >= 0.1):
                    #print "failed cut 5"
                    continue
                if( fabs(jet0_eta) >= 2.4):
                    #print "failed cut 6"
                    continue
                if (fabs(jet1_eta) >= 2.4):
                    #print "failed cut 7"
                    continue
                #if( fabs(jet2_eta) >= 2.4):
                #    #print "failed cut 6"
                #    continue
                #if (fabs(jet3_eta) >= 2.4):
                #    #print "failed cut 7"
                #    continue
                if ( corr_met <= 20):
                    #print "failed cut 8"
                    continue

                cnt += 1
                           
                if (btag0==1): nbtags += 1
                if (btag1==1): nbtags += 1

                #if (nbtags == 0): zerobtags += ((lumitype*xsectype)/neventstype)
                #if (nbtags == 1): onebtags += ((lumitype*xsectype)/neventstype)
                #if (nbtags >= 1): ge1btags += ((lumitype*xsectype)/neventstype)
                #if (nbtags >= 2): ge2btags += ((lumitype*xsectype)/neventstype)
                if (nbtags == 0): zerobtags += 1
                if (nbtags == 1): onebtags += 1
                if (nbtags >= 1): ge1btags += 1
                if (nbtags >= 2): ge2btags += 1

                syncfile.write("'"+str(run)+":"+str(lumisec)+":"+str(event)+"'\n")
                if (Jet1Jet2W_M > 2500 and nbtags >=1):
                    #highmassfile.write("'"+str(run)+":"+str(lumisec)+":"+str(event)+"', Jet1Jet2W_M = "+str(Jet1Jet2W_M)+" Jet0 pt = "+str(jet0_pt)+" Jet1 pt = "+str(jet1_pt)+" electron_pt "+str(electron_pt)+" corr_met = "+str(corr_met)+" met = "+str(met)+" BestTopMass = "+str(BestTop_M)+" \n")
                    highmassfile.write("'"+str(run)+":"+str(lumisec)+":"+str(event)+"', Jet1Jet2W_M = "+str(Jet1Jet2W_M)+" Jet0 pt = "+str(jet0_pt)+" Jet1 pt = "+str(jet1_pt)+" muon pt "+str(muon_pt)+" BestTopMass = "+str(BestTop_M)+" \n")

            print Type + ' has ' + str(Nevent) + ' total events before selection'
            print str(cnt)+' pass the cuts'            
            print str(zerobtags) + " 0 tags"
            print str(onebtags) + " 1 tags"
            print str(ge1btags) + " 1 or more tags"
            print str(ge2btags) + " 2 or more tags"


#Synch()

