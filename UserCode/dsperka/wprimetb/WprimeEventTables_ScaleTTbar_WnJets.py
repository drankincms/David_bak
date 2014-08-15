import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack, TFormula
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad

TFormula.SetMaxima(5000,5000,5000)

from array import array

from LoadData_LPC import *

global txtfile
txtfile = open('eventtables_scaleGenTopPt_WnJets.txt', 'w') 

CutList = ['preselection',
           'zerobtags',
           'onebtags',
           'ge1btags',
           'ge2btags',
           'final'
    ]

#CutList = ['preselection',
#    ]

def FillTables(channel, varName, bin, low, high):

    VariablesPre = {}
    VariablesZero = {}
    VariablesOne = {}
    VariablesGEOne = {}
    VariablesGETwo = {}
    VariablesFinal = {}
 
    WbbHist1 = {}
    WccHist1 = {}
    WbbHist2 = {}
    WccHist2 = {}
    WbbHist3 = {}
    WccHist3 = {}  
    WbbHist4 = {}
    WccHist4 = {}
 
    EventCountPre = {}
    EventCountZero = {}
    EventCountOne = {}
    EventCountGEOne = {}
    EventCountGETwo = {}
    EventCountFinal = {}
    backgroundPre = 0
    backgroundZero = 0
    backgroundOne = 0
    backgroundGEOne = 0
    backgroundGETwo = 0
    backgroundFinal = 0
   
    doqcd = 'False'

    if (channel == 'electron'):
        if doqcd == 'True':
            List = ['Data_el','ZJets_M50','WW', 'W1Jets', 'W2Jets','W3Jets','W4Jets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph','QCD_Pt_80_170_EM','QCD_Pt_170_250_EM','QCD_Pt_250_350_EM','QCD_Pt_350_EM']
        else:
            List = ['Data_el','ZJets_M50','WW', 'W1Jets', 'W2Jets','W3Jets','W4Jets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']
    if (channel == 'muon'):
        List = ['Data_mu','ZJets_M50','WW', 'W1Jets', 'W2Jets','W3Jets','W4Jets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']

    SignalList = ['Wprime800Right','Wprime900Right','Wprime1000Right','Wprime1100Right','Wprime1200Right','Wprime1300Right','Wprime1400Right','Wprime1500Right','Wprime1600Right','Wprime1700Right','Wprime1800Right','Wprime1900Right','Wprime2000Right','Wprime2100Right','Wprime2200Right','Wprime2300Right','Wprime2400Right','Wprime2500Right','Wprime2600Right','Wprime2700Right','Wprime2800Right','Wprime2900Right','Wprime3000Right']
    List.extend(SignalList)


    for cutlabel in CutList:


        if (channel == 'electron'):      
            #cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 32 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20' 
            #cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 32 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > 0.3' 
            cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 50 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20 &&  Muon_DeltaR_LjetsTopoCalcNew > 0.3' 
        if (channel == 'muon'):      
            #cut = 'jet_0_pt_WprimeCalc >= 100 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 26 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20'
            cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 50 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > 0.3'
        

        cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
        cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
        cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light

        if cutlabel == 'zerobtags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==0) '
        if cutlabel == 'onebtags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==1) '
        if cutlabel == 'ge1btags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)>=1) '
        if cutlabel == 'ge2btags': cut = cut + ' && ( (jet_0_tag_WprimeCalc==1 && (jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) || (jet_1_tag_WprimeCalc==1 && (jet_0_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) ) '
        if cutlabel == 'final':    cut = cut + ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) >= 1 ) && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 85  && Jet1Jet2_Pt_LjetsTopoCalcNew > 140 '


        #SFWjmu = 0.85             ## oldHF
        #SFWcmu = 0.92*1.66        ## oldHF
        #SFWbmu = 0.92*1.21        ## oldHF
        #SFWjmu = 1.13*0.85
        #SFWcmu = 1.24*0.92*1.66
        #SFWbmu = 1.24*0.92*1.21
        #SFWjmu = 1.05*0.85       ## myHF
        #SFWcmu = 1.25*0.92*1.66  ## myHF
        #SFWbmu = 1.25*0.92*1.21  ## myHF
        #SFWjmu = 1.08*0.85      ## myHF120
        #SFWcmu = 1.06*0.92*1.66  ## myHF120
        #SFWbmu = 1.06*0.92*1.21  ## myHF120
        #SFWjmu = 0.86*1.0   ## myHF120
        #SFWcmu = 0.95*1.66  ## myHF120
        #SFWbmu = 0.95*1.21  ## myHF120
        #SFWjmu = 0.82        ## myHF120lep50
        #SFWcmu = 0.93*1.66   ## myHF120lep50 
        #SFWbmu = 0.93*1.21   ## myHF120lep50
        #SFWjmu = 0.82        ## myHF120lep50, data driven top pt
        #SFWcmu = 0.98*1.66   ## myHF120lep50, data driven top pt
        #SFWbmu = 0.98*1.21   ## myHF120lep50, data driven top pt
        SFWjmu = 0.82        ## myHF120lep50, gen level top pt
        SFWcmu = 1.0*1.66   ## myHF120lep50, gen level top pt
        SFWbmu = 1.0*1.21   ## myHF120lep50, gen level top pt                       
        #SFWjmu = 1.0   ## noHF
        #SFWcmu = 1.66  ## noHF
        #SFWbmu = 1.21  ## noHF
      
        
        #weight_ttbar_avg = '(0.5*((1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew))+(0.9822-0.00009837*(BestJetJet2W_M_LjetsTopoCalcNew))))'
        #weight_ttbar_avg = '(1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew))'## e+jets mwp 400tp1000
        #weight_ttbar_avg = '(1.098-0.0008829*(BestTop_Pt_LjetsTopoCalcNew))'## e+jets mwp 400to750
        #weight_ttbar_avg = '(0.5*((1.098-0.0008829*(BestTop_Pt_LjetsTopoCalcNew))+(1.3-0.001659*(BestTop_Pt_LjetsTopoCalcNew))))'## e/mu+jets mwp 400to750
        #weight_ttbar_avg = '(5.809*TMath::Landau(BestTop_Pt_LjetsTopoCalcNew,149.1,154))'
        weight_ttbar_avg = 'weight_TopPt_WprimeCalc'        
                        
        if (channel=='electron'):
            #weight = 'weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc'
            weight = '( ((0.973*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
        if (channel=='muon'):
            weight = 'weight_PU_ABCD_PileUpCalc*weight_MuonEff_WprimeCalc'

        j = 0
        for Type in List:

            histName = varName+Type+cutlabel+channel
          
            WccHist1[cutlabel] = TH1D('WccHist1'+histName, 'WccHist1'+histName, bin, low, high)
            WbbHist1[cutlabel] = TH1D('WbbHist1'+histName, 'WbbHist1'+histName, bin, low, high)
            WccHist2[cutlabel] = TH1D('WccHist2'+histName, 'WccHist2'+histName, bin, low, high)
            WbbHist2[cutlabel] = TH1D('WbbHist2'+histName, 'WbbHist2'+histName, bin, low, high)
            WccHist3[cutlabel] = TH1D('WccHist3'+histName, 'WccHist3'+histName, bin, low, high)
            WbbHist3[cutlabel] = TH1D('WbbHist3'+histName, 'WbbHist3'+histName, bin, low, high)
            WccHist4[cutlabel] = TH1D('WccHist4'+histName, 'WccHist4'+histName, bin, low, high)
            WbbHist4[cutlabel] = TH1D('WbbHist4'+histName, 'WbbHist4'+histName, bin, low, high)

            SF = 1.0    

            if cutlabel == 'preselection': VariablesPre[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'zerobtags': VariablesZero[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'onebtags': VariablesOne[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'ge1btags': VariablesGEOne[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'ge2btags': VariablesGETwo[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'final': VariablesFinal[Type] = TH1D(histName, histName, bin, low, high)
               
            if Type.startswith('Data'):
                Trees[Type].Draw(var + " >> " + histName, "(" + cut + ")", 'goff')
            elif Type == 'W1Jets':       
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WbbHist1"+histName, "("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WccHist1"+histName, "("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
            elif Type == 'W2Jets':       
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WbbHist2"+histName, "("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WccHist2"+histName, "("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
            elif Type == 'W3Jets':       
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WbbHist3"+histName, "("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WccHist3"+histName, "("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
            elif Type == 'W4Jets':       
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WbbHist4"+histName, "("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WccHist4"+histName, "("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
                
            elif Type.startswith('TTbar'):
                Trees[Type].Draw(var + " >> " + histName, "("+weight_ttbar_avg+")*("+weight+")*(" + cut + ")", 'goff')
            else:
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + cut + ")", 'goff')

            if (channel == 'electron'):
                lumi = lumi_el
            if (channel == 'muon'):
                lumi = lumi_mu
            
            if (not Type.startswith('Data')):
                if Type.startswith('W1Jets'): 
                    WbbHist1[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    WccHist1[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if Type.startswith('W2Jets'): 
                    WbbHist2[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    WccHist2[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if Type.startswith('W3Jets'): 
                    WbbHist3[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    WccHist3[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if Type.startswith('W4Jets'): 
                    WbbHist4[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    WccHist4[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'preselection':
                    if VariablesPre[Type].Integral() != 0:
                       VariablesPre[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'zerobtags':
                    if VariablesZero[Type].Integral() != 0:
                       VariablesZero[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'onebtags':
                    if VariablesOne[Type].Integral() != 0:
                       VariablesOne[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'ge1btags':
                    if VariablesGEOne[Type].Integral() != 0:
                       VariablesGEOne[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'ge2btags':
                    if VariablesGETwo[Type].Integral() != 0:
                       VariablesGETwo[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'final':
                    if VariablesFinal[Type].Integral() != 0:
                       VariablesFinal[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'preselection':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundPre += VariablesPre[Type].Integral()
                    EventCountPre[Type] = VariablesPre[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountPre['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountPre['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundPre += WbbHist1[cutlabel].Integral()
                        backgroundPre += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountPre['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountPre['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundPre += WbbHist2[cutlabel].Integral()
                        backgroundPre += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountPre['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountPre['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundPre += WbbHist3[cutlabel].Integral()
                        backgroundPre += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountPre['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountPre['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundPre += WbbHist4[cutlabel].Integral()
                        backgroundPre += WccHist4[cutlabel].Integral()
                if cutlabel == 'zerobtags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundZero += VariablesZero[Type].Integral()
                    EventCountZero[Type] = VariablesZero[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountZero['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountZero['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundZero += WbbHist1[cutlabel].Integral()
                        backgroundZero += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountZero['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountZero['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundZero += WbbHist2[cutlabel].Integral()
                        backgroundZero += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountZero['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountZero['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundZero += WbbHist3[cutlabel].Integral()
                        backgroundZero += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountZero['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountZero['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundZero += WbbHist4[cutlabel].Integral()
                        backgroundZero += WccHist4[cutlabel].Integral()
                if cutlabel == 'onebtags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundOne += VariablesOne[Type].Integral()
                    EventCountOne[Type] = VariablesOne[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountOne['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountOne['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundOne += WbbHist1[cutlabel].Integral()
                        backgroundOne += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountOne['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountOne['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundOne += WbbHist2[cutlabel].Integral()
                        backgroundOne += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountOne['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountOne['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundOne += WbbHist3[cutlabel].Integral()
                        backgroundOne += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountOne['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountOne['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundOne += WbbHist4[cutlabel].Integral()
                        backgroundOne += WccHist4[cutlabel].Integral()
                if cutlabel == 'ge1btags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundGEOne += VariablesGEOne[Type].Integral()
                    EventCountGEOne[Type] = VariablesGEOne[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountGEOne['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountGEOne['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundGEOne += WbbHist1[cutlabel].Integral()
                        backgroundGEOne += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountGEOne['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountGEOne['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundGEOne += WbbHist2[cutlabel].Integral()
                        backgroundGEOne += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountGEOne['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountGEOne['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundGEOne += WbbHist3[cutlabel].Integral()
                        backgroundGEOne += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountGEOne['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountGEOne['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundGEOne += WbbHist4[cutlabel].Integral()
                        backgroundGEOne += WccHist4[cutlabel].Integral()
                if cutlabel == 'ge2btags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundGETwo += VariablesGETwo[Type].Integral()
                    EventCountGETwo[Type] = VariablesGETwo[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountGETwo['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountGETwo['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundGETwo += WbbHist1[cutlabel].Integral()
                        backgroundGETwo += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountGETwo['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountGETwo['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundGETwo += WbbHist2[cutlabel].Integral()
                        backgroundGETwo += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountGETwo['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountGETwo['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundGETwo += WbbHist3[cutlabel].Integral()
                        backgroundGETwo += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountGETwo['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountGETwo['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundGETwo += WbbHist4[cutlabel].Integral()
                        backgroundGETwo += WccHist4[cutlabel].Integral()
                if cutlabel == 'final':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundFinal += VariablesFinal[Type].Integral()
                    EventCountFinal[Type] = VariablesFinal[Type].Integral()
                    if Type.startswith('W1Jets'): 
                        EventCountFinal['Wbb1'] = WbbHist1[cutlabel].Integral()
                        EventCountFinal['Wcc1'] = WccHist1[cutlabel].Integral()
                        backgroundFinal += WbbHist1[cutlabel].Integral()
                        backgroundFinal += WccHist1[cutlabel].Integral()
                    if Type.startswith('W2Jets'): 
                        EventCountFinal['Wbb2'] = WbbHist2[cutlabel].Integral()
                        EventCountFinal['Wcc2'] = WccHist2[cutlabel].Integral()
                        backgroundFinal += WbbHist2[cutlabel].Integral()
                        backgroundFinal += WccHist2[cutlabel].Integral()
                    if Type.startswith('W3Jets'): 
                        EventCountFinal['Wbb3'] = WbbHist3[cutlabel].Integral()
                        EventCountFinal['Wcc3'] = WccHist3[cutlabel].Integral()
                        backgroundFinal += WbbHist3[cutlabel].Integral()
                        backgroundFinal += WccHist3[cutlabel].Integral()
                    if Type.startswith('W4Jets'): 
                        EventCountFinal['Wbb4'] = WbbHist4[cutlabel].Integral()
                        EventCountFinal['Wcc4'] = WccHist4[cutlabel].Integral()
                        backgroundFinal += WbbHist4[cutlabel].Integral()
                        backgroundFinal += WccHist4[cutlabel].Integral()
                    
            if Type.startswith('Data'):
                if cutlabel == 'preselection':
                    EventCountPre[Type] = VariablesPre[Type].Integral()
                if cutlabel == 'zerobtags':
                    EventCountZero[Type] = VariablesZero[Type].Integral()
                if cutlabel == 'onebtags':
                    EventCountOne[Type] = VariablesOne[Type].Integral()
                if cutlabel == 'ge1btags':
                    EventCountGEOne[Type] = VariablesGEOne[Type].Integral()
                if cutlabel == 'ge2btags':
                    EventCountGETwo[Type] = VariablesGETwo[Type].Integral()
                if cutlabel == 'final':
                    EventCountFinal[Type] = VariablesFinal[Type].Integral()
            j=j+1


    print channel,' channel >=1 b-tag W+jets counts:'
    print 'W+1Jets ',str(int(EventCountGEOne['W1Jets']+EventCountGEOne['Wbb1']+EventCountGEOne['Wcc1']))
    print 'W+2Jets ',str(int(EventCountGEOne['W2Jets']+EventCountGEOne['Wbb2']+EventCountGEOne['Wcc2']))
    print 'W+3Jets ',str(int(EventCountGEOne['W3Jets']+EventCountGEOne['Wbb3']+EventCountGEOne['Wcc3']))
    print 'W+4Jets ',str(int(EventCountGEOne['W4Jets']+EventCountGEOne['Wbb4']+EventCountGEOne['Wcc4']))
                             
    print channel,' channel ==1 b-tag W+jets counts:'
    print 'W+1Jets ',str(int(EventCountOne['W1Jets']+EventCountOne['Wbb1']+EventCountOne['Wcc1']))
    print 'W+2Jets ',str(int(EventCountOne['W2Jets']+EventCountOne['Wbb2']+EventCountOne['Wcc2']))
    print 'W+3Jets ',str(int(EventCountOne['W3Jets']+EventCountOne['Wbb3']+EventCountOne['Wcc3']))
    print 'W+4Jets ',str(int(EventCountOne['W4Jets']+EventCountOne['Wbb4']+EventCountOne['Wcc4']))


                    
    EventCountPre['WJets'] = EventCountPre['W1Jets']+EventCountPre['W2Jets']+EventCountPre['W3Jets']+EventCountPre['W4Jets']
    EventCountZero['WJets'] = EventCountZero['W1Jets']+EventCountZero['W2Jets']+EventCountZero['W3Jets']+EventCountZero['W4Jets']
    EventCountOne['WJets'] = EventCountOne['W1Jets']+EventCountOne['W2Jets']+EventCountOne['W3Jets']+EventCountOne['W4Jets']
    EventCountGEOne['WJets'] = EventCountGEOne['W1Jets']+EventCountGEOne['W2Jets']+EventCountGEOne['W3Jets']+EventCountGEOne['W4Jets']
    EventCountGETwo['WJets'] = EventCountGETwo['W1Jets']+EventCountGETwo['W2Jets']+EventCountGETwo['W3Jets']+EventCountGETwo['W4Jets']
    EventCountFinal['WJets'] = EventCountFinal['W1Jets']+EventCountFinal['W2Jets']+EventCountFinal['W3Jets']+EventCountFinal['W4Jets']

    EventCountPre['Wbb'] = EventCountPre['Wbb1']+EventCountPre['Wbb2']+EventCountPre['Wbb3']+EventCountPre['Wbb4']
    EventCountZero['Wbb'] = EventCountZero['Wbb1']+EventCountZero['Wbb2']+EventCountZero['Wbb3']+EventCountZero['Wbb4']
    EventCountOne['Wbb'] = EventCountOne['Wbb1']+EventCountOne['Wbb2']+EventCountOne['Wbb3']+EventCountOne['Wbb4']
    EventCountGEOne['Wbb'] = EventCountGEOne['Wbb1']+EventCountGEOne['Wbb2']+EventCountGEOne['Wbb3']+EventCountGEOne['Wbb4']
    EventCountGETwo['Wbb'] = EventCountGETwo['Wbb1']+EventCountGETwo['Wbb2']+EventCountGETwo['Wbb3']+EventCountGETwo['Wbb4']
    EventCountFinal['Wbb'] = EventCountFinal['Wbb1']+EventCountFinal['Wbb2']+EventCountFinal['Wbb3']+EventCountFinal['Wbb4']

    EventCountPre['Wcc'] = EventCountPre['Wcc1']+EventCountPre['Wcc2']+EventCountPre['Wcc3']+EventCountPre['Wcc4']
    EventCountZero['Wcc'] = EventCountZero['Wcc1']+EventCountZero['Wcc2']+EventCountZero['Wcc3']+EventCountZero['Wcc4']
    EventCountOne['Wcc'] = EventCountOne['Wcc1']+EventCountOne['Wcc2']+EventCountOne['Wcc3']+EventCountOne['Wcc4']
    EventCountGEOne['Wcc'] = EventCountGEOne['Wcc1']+EventCountGEOne['Wcc2']+EventCountGEOne['Wcc3']+EventCountGEOne['Wcc4']
    EventCountGETwo['Wcc'] = EventCountGETwo['Wcc1']+EventCountGETwo['Wcc2']+EventCountGETwo['Wcc3']+EventCountGETwo['Wcc4']
    EventCountFinal['Wcc'] = EventCountFinal['Wcc1']+EventCountFinal['Wcc2']+EventCountFinal['Wcc3']+EventCountFinal['Wcc4']

    if (channel == 'electron' and doqcd == 'True'):
        EventCountPre['QCD_Pt_80_170_EM'] += EventCountPre['QCD_Pt_170_250_EM']
        EventCountPre['QCD_Pt_80_170_EM'] += EventCountPre['QCD_Pt_250_350_EM']
        EventCountPre['QCD_Pt_80_170_EM'] += EventCountPre['QCD_Pt_350_EM']
      
        EventCountZero['QCD_Pt_80_170_EM'] += EventCountZero['QCD_Pt_170_250_EM']
        EventCountZero['QCD_Pt_80_170_EM'] += EventCountZero['QCD_Pt_250_350_EM']
        EventCountZero['QCD_Pt_80_170_EM'] += EventCountZero['QCD_Pt_350_EM']
      
        EventCountOne['QCD_Pt_80_170_EM'] += EventCountOne['QCD_Pt_170_250_EM']
        EventCountOne['QCD_Pt_80_170_EM'] += EventCountOne['QCD_Pt_250_350_EM']
        EventCountOne['QCD_Pt_80_170_EM'] += EventCountOne['QCD_Pt_350_EM']
      
        EventCountGEOne['QCD_Pt_80_170_EM'] += EventCountGEOne['QCD_Pt_170_250_EM']
        EventCountGEOne['QCD_Pt_80_170_EM'] += EventCountGEOne['QCD_Pt_250_350_EM']
        EventCountGEOne['QCD_Pt_80_170_EM'] += EventCountGEOne['QCD_Pt_350_EM']
      
        EventCountGETwo['QCD_Pt_80_170_EM'] += EventCountGETwo['QCD_Pt_170_250_EM']
        EventCountGETwo['QCD_Pt_80_170_EM'] += EventCountGETwo['QCD_Pt_250_350_EM']
        EventCountGETwo['QCD_Pt_80_170_EM'] += EventCountGETwo['QCD_Pt_350_EM']
      
        EventCountFinal['QCD_Pt_80_170_EM'] += EventCountFinal['QCD_Pt_170_250_EM']
        EventCountFinal['QCD_Pt_80_170_EM'] += EventCountFinal['QCD_Pt_250_350_EM']
        EventCountFinal['QCD_Pt_80_170_EM'] += EventCountFinal['QCD_Pt_350_EM']
      

    if (channel == 'electron'): suffix = 'el'
    if (channel == 'muon'): suffix = 'mu'

    txtfile.write("\\begin{table}[!h!tb]\n")
    txtfile.write("\\begin{center}\n")
    txtfile.write("\small\n")
    txtfile.write("      \caption{\n")
    txtfile.write("	Number of selected data, and background  events in the "+channel+" channel. \n")
    txtfile.write("For the background samples, the expectation is computed\n")
    txtfile.write("corresponding to an integrated luminosity of $"+str(lumi)+" \pbinv$.\n")
    txtfile.write("\label{tab:"+suffix+"_cut_flow_bkg}\n") 
    txtfile.write("        }\n")
    txtfile.write("\\begin{tabular}{l|c|c|c|c|c|c}\n")
    txtfile.write("\hline\n")
    txtfile.write("Process    & \multicolumn{6}{c}{Number of Events} \\\ \hline\hline\n")
    txtfile.write("   & & \multicolumn{4}{c|}{Number of b-tagged jets} & $p_{T}^{top}>$ 85,$p_{T}^{jet1,jet2}>$140,130$<m_{top}<$210 \\\ \n")
    txtfile.write(" & Pre selection & = 0 & = 1 \n")
    txtfile.write("& $>$ 0 & $>$ 1  & $>$0 b-tagged jet    \\\ \hline\n")
    txtfile.write("{\\bf Data} & "+str(int(round(EventCountPre['Data_'+suffix])))+" & "+str(int(round(EventCountZero['Data_'+suffix])))+" & "+str(int(round(EventCountOne['Data_'+suffix])))+" & "+str(int(round(EventCountGEOne['Data_'+suffix])))+" & "+str(int(round(EventCountGETwo['Data_'+suffix])))+" & "+str(int(round(EventCountFinal['Data_'+suffix])))+"   \\\ \hline\n") 
    txtfile.write("\multicolumn{6}{l}{\\bf Background:}     \\\ \hline\n")
    txtfile.write("$t\overline t$ &"+str(int(round(EventCountPre['TTbar_Madgraph'])))+" & "+str(int(round(EventCountZero['TTbar_Madgraph'])))+" & "+str(int(round(EventCountOne['TTbar_Madgraph'])))+" & "+str(int(round(EventCountGEOne['TTbar_Madgraph'])))+" & "+str(int(round(EventCountGETwo['TTbar_Madgraph'])))+" & "+str(int(round(EventCountFinal['TTbar_Madgraph'])))+"  \\\ \n")
    #txtfile.write("$t\overline t$ &"+str(int(round(EventCountPre['TTbar_Powheg'])))+" & "+str(int(round(EventCountZero['TTbar_Powheg'])))+" & "+str(int(round(EventCountOne['TTbar_Powheg'])))+" & "+str(int(round(EventCountGEOne['TTbar_Powheg'])))+" & "+str(int(round(EventCountGETwo['TTbar_Powheg'])))+" & "+str(int(round(EventCountFinal['TTbar_Powheg'])))+"  \\\ \n")
    txtfile.write("$tqb$ &"+str(int(round(EventCountPre['T_t'])))+" & "+str(int(round(EventCountZero['T_t'])))+" & "+str(int(round(EventCountOne['T_t'])))+" & "+str(int(round(EventCountGEOne['T_t'])))+" & "+str(int(round(EventCountGETwo['T_t'])))+"  & "+str(int(round(EventCountFinal['T_t'])))+"  \\\ \n")
    txtfile.write("$\overline{t}qb$ &"+str(int(round(EventCountPre['Tbar_t'])))+" & "+str(int(round(EventCountZero['Tbar_t'])))+" & "+str(int(round(EventCountOne['Tbar_t'])))+" & "+str(int(round(EventCountGEOne['Tbar_t'])))+" & "+str(int(round(EventCountGETwo['Tbar_t'])))+"  & "+str(int(round(EventCountFinal['Tbar_t'])))+"  \\\ \n")
    txtfile.write("$tW$  &"+str(int(round(EventCountPre['T_tW'])))+" & "+str(int(round(EventCountZero['T_tW'])))+" & "+str(int(round(EventCountOne['T_tW'])))+" & "+str(int(round(EventCountGEOne['T_tW'])))+" & "+str(int(round(EventCountGETwo['T_tW'])))+"  & "+str(int(round(EventCountFinal['T_tW'])))+"  \\\ \n")
    txtfile.write("$\overline{t}W$ &"+str(int(round(EventCountPre['Tbar_tW'])))+" & "+str(int(round(EventCountZero['Tbar_tW'])))+" & "+str(int(round(EventCountOne['Tbar_tW'])))+" & "+str(int(round(EventCountGEOne['Tbar_tW'])))+" & "+str(int(round(EventCountGETwo['Tbar_tW'])))+"  & "+str(int(round(EventCountFinal['Tbar_tW'])))+"  \\\ \n")
    txtfile.write("$tb$ & "+str(int(round(EventCountPre['T_s'])))+" & "+str(int(round(EventCountZero['T_s'])))+" & "+str(int(round(EventCountOne['T_s'])))+" & "+str(int(round(EventCountGEOne['T_s'])))+" & "+str(int(round(EventCountGETwo['T_s'])))+" & "+str(int(round(EventCountFinal['T_s'])))+"    \\\ \n")
    txtfile.write("$\overline{t}b$ & "+str(int(round(EventCountPre['Tbar_s'])))+" & "+str(int(round(EventCountZero['Tbar_s'])))+" & "+str(int(round(EventCountOne['Tbar_s'])))+" & "+str(int(round(EventCountGEOne['Tbar_s'])))+" & "+str(int(round(EventCountGETwo['Tbar_s'])))+" & "+str(int(round(EventCountFinal['Tbar_s'])))+"    \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+light jets &"+str(int(round(EventCountPre['WJets'])))+" & "+str(int(round(EventCountZero['WJets'])))+" & "+str(int(round(EventCountOne['WJets'])))+" & "+str(int(round(EventCountGEOne['WJets'])))+" & "+str(int(round(EventCountGETwo['WJets'])))+"  & "+str(int(round(EventCountFinal['WJets'])))+"  \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+c jets &"+str(int(round(EventCountPre['Wcc'])))+" & "+str(int(round(EventCountZero['Wcc'])))+" & "+str(int(round(EventCountOne['Wcc'])))+" & "+str(int(round(EventCountGEOne['Wcc'])))+" & "+str(int(round(EventCountGETwo['Wcc'])))+"  & "+str(int(round(EventCountFinal['Wcc'])))+"  \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+b jets &"+str(int(round(EventCountPre['Wbb'])))+" & "+str(int(round(EventCountZero['Wbb'])))+" & "+str(int(round(EventCountOne['Wbb'])))+" & "+str(int(round(EventCountGEOne['Wbb'])))+" & "+str(int(round(EventCountGETwo['Wbb'])))+"  & "+str(int(round(EventCountFinal['Wbb'])))+"  \\\ \n")
    txtfile.write("$Z/\gamma^*(\\rightarrow\ell\ell$)+jets &"+str(int(round(EventCountPre['ZJets_M50'])))+" & "+str(int(round(EventCountZero['ZJets_M50'])))+" & "+str(int(round(EventCountOne['ZJets_M50'])))+" & "+str(int(round(EventCountGEOne['ZJets_M50'])))+" & "+str(int(round(EventCountGETwo['ZJets_M50'])))+"  & "+str(int(round(EventCountFinal['ZJets_M50'])))+" \\\ \n")
    txtfile.write("$WW$ &"+str(int(round(EventCountPre['WW'])))+" & "+str(int(round(EventCountZero['WW'])))+" & "+str(int(round(EventCountOne['WW'])))+" & "+str(int(round(EventCountGEOne['WW'])))+" & "+str(int(round(EventCountGETwo['WW'])))+"  & "+str(int(round(EventCountFinal['WW'])))+" \\\ \n")
    if (channel=='electron' and doqcd == 'True'):
        txtfile.write("$QCD$ &"+str(int(round(EventCountPre['QCD_Pt_80_170_EM'])))+" & "+str(int(round(EventCountZero['QCD_Pt_80_170_EM'])))+" & "+str(int(round(EventCountOne['QCD_Pt_80_170_EM'])))+" & "+str(int(round(EventCountGEOne['QCD_Pt_80_170_EM'])))+" & "+str(int(round(EventCountGETwo['QCD_Pt_80_170_EM'])))+"  & "+str(int(round(EventCountFinal['QCD_Pt_80_170_EM'])))+" \\\ \n")
    txtfile.write("{\\bf Total Background} &"+str(int(round(backgroundPre)))+" & "+str(int(round(backgroundZero)))+" & "+str(int(round(backgroundOne)))+" & "+str(int(round(backgroundGEOne)))+" & "+str(int(round(backgroundGETwo)))+" & "+str(int(round(backgroundFinal)))+" \\\ \hline\n")
    txtfile.write("{\\bf MC / Data} &"+str(float(round(backgroundPre/EventCountPre['Data_'+suffix], 3)))+" & "+str(float(round(backgroundZero/EventCountZero['Data_'+suffix], 3)))+" & "+str(float(round(backgroundOne/EventCountOne['Data_'+suffix], 3)))+" & "+str(float(round(backgroundGEOne/EventCountGEOne['Data_'+suffix], 3)))+" & "+str(float(round(backgroundGETwo/EventCountGETwo['Data_'+suffix], 3)))+" & "+str(float(round(backgroundFinal/EventCountFinal['Data_'+suffix], 3)))+" \\\ \hline\n")
    txtfile.write("\hline\n")
    txtfile.write("\end{tabular}\n")
    txtfile.write("\\normalsize\n")
    txtfile.write("\end{center}\n")
    txtfile.write("\end{table}\n")
    txtfile.write("\n")
    txtfile.write("\n")
    txtfile.write("\n")

    EventCountPre['T_s'] += EventCountPre['Tbar_s']      
    EventCountZero['T_s'] += EventCountZero['Tbar_s'] 
    EventCountOne['T_s'] += EventCountOne['Tbar_s']  
    EventCountGEOne['T_s'] += EventCountGEOne['Tbar_s'] 
    EventCountGETwo['T_s'] += EventCountGETwo['Tbar_s'] 
    EventCountFinal['T_s'] += EventCountFinal['Tbar_s'] 
      
    txtfile.write("\\begin{table}[!h!tb]\n")
    txtfile.write("\\begin{center}\n")
    txtfile.write("\small\n")
    txtfile.write("     \caption{\n")
    txtfile.write("       Number of selected  signal events in the "+channel+" channel.\n")
    txtfile.write("For the  signal samples, the expectation is computed\n")
    txtfile.write("corresponding to an integrated luminosity of $"+str(lumi)+" \pbinv$.\n")
    txtfile.write("\label{tab:"+suffix+"_cut_flow_Wprime}\n")
    txtfile.write("        }\n")
    txtfile.write("\\begin{tabular}{l|c|c|c|c|c|c} \n")
    txtfile.write("\hline \n")
    txtfile.write("Process    & \multicolumn{6}{c}{Number of Events} \\\ \hline\hline \n")
    txtfile.write("   & & \multicolumn{4}{c|}{Number of b-tagged jets} & $p_{T}^{top}>$ 85,$p_{T}^{jet1,jet2}>$140 \\\ \n")
    txtfile.write(" & Pre selection & = 0 & = 1\n")
    txtfile.write("& $>$ 0 & $>$ 1 & $>$0 b-tagged jet      \\\ \hline \n")
    txtfile.write("\multicolumn{6}{l}{\\bf Signal:}     \\\ \hline \n")
    txtfile.write("$tb$ & "+str(int(round(EventCountPre['T_s'])))+" & "+str(int(round(EventCountZero['T_s'])))+" & "+str(int(round(EventCountOne['T_s'])))+" & "+str(int(round(EventCountGEOne['T_s'])))+" & "+str(int(round(EventCountGETwo['T_s'])))+" & "+str(int(round(EventCountFinal['T_s'])))+"    \\\ \n")
    txtfile.write("{\underline {\\bf ${\PWpr}_R\\rightarrow tb$}} & &  &  & & \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 800 GeV &"+str(int(round(EventCountPre['Wprime800Right'])))+" & "+str(int(round(EventCountZero['Wprime800Right'])))+" & "+str(int(round(EventCountOne['Wprime800Right'])))+" & "+str(int(round(EventCountGEOne['Wprime800Right'])))+" & "+str(int(round(EventCountGETwo['Wprime800Right'])))+" & "+str(int(round(EventCountFinal['Wprime800Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 900 GeV &"+str(int(round(EventCountPre['Wprime900Right'])))+" & "+str(int(round(EventCountZero['Wprime900Right'])))+" & "+str(int(round(EventCountOne['Wprime900Right'])))+" & "+str(int(round(EventCountGEOne['Wprime900Right'])))+" & "+str(int(round(EventCountGETwo['Wprime900Right'])))+" & "+str(int(round(EventCountFinal['Wprime900Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1000 GeV &"+str(int(round(EventCountPre['Wprime1000Right'])))+" & "+str(int(round(EventCountZero['Wprime1000Right'])))+" & "+str(int(round(EventCountOne['Wprime1000Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1000Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1000Right'])))+" & "+str(int(round(EventCountFinal['Wprime1000Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1100 GeV &"+str(int(round(EventCountPre['Wprime1100Right'])))+" & "+str(int(round(EventCountZero['Wprime1100Right'])))+" & "+str(int(round(EventCountOne['Wprime1100Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1100Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1100Right'])))+" & "+str(int(round(EventCountFinal['Wprime1100Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1200 GeV &"+str(int(round(EventCountPre['Wprime1200Right'])))+" & "+str(int(round(EventCountZero['Wprime1200Right'])))+" & "+str(int(round(EventCountOne['Wprime1200Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1200Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1200Right'])))+" & "+str(int(round(EventCountFinal['Wprime1200Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1300 GeV &"+str(int(round(EventCountPre['Wprime1300Right'])))+" & "+str(int(round(EventCountZero['Wprime1300Right'])))+" & "+str(int(round(EventCountOne['Wprime1300Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1300Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1300Right'])))+" & "+str(int(round(EventCountFinal['Wprime1300Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1400 GeV &"+str(int(round(EventCountPre['Wprime1400Right'])))+" & "+str(int(round(EventCountZero['Wprime1400Right'])))+" & "+str(int(round(EventCountOne['Wprime1400Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1400Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1400Right'])))+" & "+str(int(round(EventCountFinal['Wprime1400Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1500 GeV &"+str(int(round(EventCountPre['Wprime1500Right'])))+" & "+str(int(round(EventCountZero['Wprime1500Right'])))+" & "+str(int(round(EventCountOne['Wprime1500Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1500Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1500Right'])))+" & "+str(int(round(EventCountFinal['Wprime1500Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1600 GeV &"+str(int(round(EventCountPre['Wprime1600Right'])))+" & "+str(int(round(EventCountZero['Wprime1600Right'])))+" & "+str(int(round(EventCountOne['Wprime1600Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1600Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1600Right'])))+" & "+str(int(round(EventCountFinal['Wprime1600Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1700 GeV &"+str(int(round(EventCountPre['Wprime1700Right'])))+" & "+str(int(round(EventCountZero['Wprime1700Right'])))+" & "+str(int(round(EventCountOne['Wprime1700Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1700Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1700Right'])))+" & "+str(int(round(EventCountFinal['Wprime1700Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1800 GeV &"+str(int(round(EventCountPre['Wprime1800Right'])))+" & "+str(int(round(EventCountZero['Wprime1800Right'])))+" & "+str(int(round(EventCountOne['Wprime1800Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1800Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1800Right'])))+" & "+str(int(round(EventCountFinal['Wprime1800Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1900 GeV &"+str(int(round(EventCountPre['Wprime1900Right'])))+" & "+str(int(round(EventCountZero['Wprime1900Right'])))+" & "+str(int(round(EventCountOne['Wprime1900Right'])))+" & "+str(int(round(EventCountGEOne['Wprime1900Right'])))+" & "+str(int(round(EventCountGETwo['Wprime1900Right'])))+" & "+str(int(round(EventCountFinal['Wprime1900Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2000 GeV &"+str(int(round(EventCountPre['Wprime2000Right'])))+" & "+str(int(round(EventCountZero['Wprime2000Right'])))+" & "+str(int(round(EventCountOne['Wprime2000Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2000Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2000Right'])))+" & "+str(int(round(EventCountFinal['Wprime2000Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2100 GeV &"+str(int(round(EventCountPre['Wprime2100Right'])))+" & "+str(int(round(EventCountZero['Wprime2100Right'])))+" & "+str(int(round(EventCountOne['Wprime2100Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2100Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2100Right'])))+" & "+str(int(round(EventCountFinal['Wprime2100Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2200 GeV &"+str(int(round(EventCountPre['Wprime2200Right'])))+" & "+str(int(round(EventCountZero['Wprime2200Right'])))+" & "+str(int(round(EventCountOne['Wprime2200Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2200Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2200Right'])))+" & "+str(int(round(EventCountFinal['Wprime2200Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2300 GeV &"+str(int(round(EventCountPre['Wprime2300Right'])))+" & "+str(int(round(EventCountZero['Wprime2300Right'])))+" & "+str(int(round(EventCountOne['Wprime2300Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2300Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2300Right'])))+" & "+str(int(round(EventCountFinal['Wprime2300Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2500 GeV &"+str(int(round(EventCountPre['Wprime2500Right'])))+" & "+str(int(round(EventCountZero['Wprime2500Right'])))+" & "+str(int(round(EventCountOne['Wprime2500Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2500Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2500Right'])))+" & "+str(int(round(EventCountFinal['Wprime2500Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2400 GeV &"+str(int(round(EventCountPre['Wprime2400Right'])))+" & "+str(int(round(EventCountZero['Wprime2400Right'])))+" & "+str(int(round(EventCountOne['Wprime2400Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2400Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2400Right'])))+" & "+str(int(round(EventCountFinal['Wprime2400Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2600 GeV &"+str(int(round(EventCountPre['Wprime2600Right'])))+" & "+str(int(round(EventCountZero['Wprime2600Right'])))+" & "+str(int(round(EventCountOne['Wprime2600Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2600Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2600Right'])))+" & "+str(int(round(EventCountFinal['Wprime2600Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2700 GeV &"+str(int(round(EventCountPre['Wprime2700Right'])))+" & "+str(int(round(EventCountZero['Wprime2700Right'])))+" & "+str(int(round(EventCountOne['Wprime2700Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2700Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2700Right'])))+" & "+str(int(round(EventCountFinal['Wprime2700Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2800 GeV &"+str(int(round(EventCountPre['Wprime2800Right'])))+" & "+str(int(round(EventCountZero['Wprime2800Right'])))+" & "+str(int(round(EventCountOne['Wprime2800Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2800Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2800Right'])))+" & "+str(int(round(EventCountFinal['Wprime2800Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2900 GeV &"+str(int(round(EventCountPre['Wprime2900Right'])))+" & "+str(int(round(EventCountZero['Wprime2900Right'])))+" & "+str(int(round(EventCountOne['Wprime2900Right'])))+" & "+str(int(round(EventCountGEOne['Wprime2900Right'])))+" & "+str(int(round(EventCountGETwo['Wprime2900Right'])))+" & "+str(int(round(EventCountFinal['Wprime2900Right'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 3000 GeV &"+str(int(round(EventCountPre['Wprime3000Right'])))+" & "+str(int(round(EventCountZero['Wprime3000Right'])))+" & "+str(int(round(EventCountOne['Wprime3000Right'])))+" & "+str(int(round(EventCountGEOne['Wprime3000Right'])))+" & "+str(int(round(EventCountGETwo['Wprime3000Right'])))+" & "+str(int(round(EventCountFinal['Wprime3000Right'])))+"   \\\ \n")

    '''
    txtfile.write("{\underline {\\bf ${\PWpr}_L\\rightarrow tb$}}  &  & &  & & \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)800 GeV &"+str(int(round(EventCountPre['Wprime_800_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_800_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_800_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_800_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_800_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_800_LeftWprime'])))+"    \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1000 GeV &"+str(int(round(EventCountPre['Wprime_1000_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1000_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1000_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1000_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1000_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1000_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1100 GeV &"+str(int(round(EventCountPre['Wprime_1100_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1100_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1100_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1100_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1100_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1100_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1200 GeV &"+str(int(round(EventCountPre['Wprime_1200_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1200_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1200_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1200_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1200_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1200_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1300 GeV &"+str(int(round(EventCountPre['Wprime_1300_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1300_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1300_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1300_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1300_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1300_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1400 GeV &"+str(int(round(EventCountPre['Wprime_1400_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1400_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1400_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1400_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1400_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1400_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1500 GeV &"+str(int(round(EventCountPre['Wprime_1500_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1500_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1500_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1500_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1500_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1500_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1600 GeV &"+str(int(round(EventCountPre['Wprime_1600_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1600_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1600_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1600_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1600_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1600_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1700 GeV &"+str(int(round(EventCountPre['Wprime_1700_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1700_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1700_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1700_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1700_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1700_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1900 GeV &"+str(int(round(EventCountPre['Wprime_1900_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1900_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1900_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1900_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1900_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1900_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)2100 GeV &"+str(int(round(EventCountPre['Wprime_2100_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_2100_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2100_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2100_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2100_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2100_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)2300 GeV &"+str(int(round(EventCountPre['Wprime_2300_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_2300_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2300_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2300_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2300_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2300_LeftWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)2500 GeV &"+str(int(round(EventCountPre['Wprime_2500_LeftWprime'])))+" & "+str(int(round(EventCountZero['Wprime_2500_LeftWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2500_LeftWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2500_LeftWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2500_LeftWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2500_LeftWprime'])))+"   \\\ \n")
    txtfile.write("{\underline {\\bf ${\PWpr}_{LR}\\rightarrow tb$}}&  & & & & \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)800 GeV &"+str(int(round(EventCountPre['Wprime_800_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_800_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_800_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_800_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_800_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_800_MixRLWprime'])))+"    \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)900 GeV &"+str(int(round(EventCountPre['Wprime_900_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_900_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_900_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_900_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_900_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_900_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1000 GeV &"+str(int(round(EventCountPre['Wprime_1000_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1000_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1000_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1000_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1000_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1000_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1100 GeV &"+str(int(round(EventCountPre['Wprime_1100_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1100_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1100_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1100_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1100_MixRLWprime'])))+"  & "+str(int(round(EventCountFinal['Wprime_1100_MixRLWprime'])))+"  \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1200 GeV & "+str(int(round(EventCountPre['Wprime_1200_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1200_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1200_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1200_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1200_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1200_MixRLWprime'])))+"  \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1300 GeV &"+str(int(round(EventCountPre['Wprime_1300_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1300_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1300_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1300_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1300_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1300_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1400 GeV &"+str(int(round(EventCountPre['Wprime_1400_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1400_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1400_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1400_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1400_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1400_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1500 GeV &"+str(int(round(EventCountPre['Wprime_1500_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1500_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1500_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1500_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1500_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1500_MixRLWprime'])))+"    \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1600 GeV &"+str(int(round(EventCountPre['Wprime_1600_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1600_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1600_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1600_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1600_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1600_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1700 GeV &"+str(int(round(EventCountPre['Wprime_1700_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1700_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1700_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1700_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1700_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1700_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1900 GeV &"+str(int(round(EventCountPre['Wprime_1900_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_1900_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_1900_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_1900_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_1900_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_1900_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2100 GeV &"+str(int(round(EventCountPre['Wprime_2100_MixRLWprime'])))+" & "+str(int(round(EventCountZero['Wprime_2100_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2100_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2100_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2100_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2100_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2300 GeV &"+str(EventCountPre['Wprime_2300_MixRLWprime'])+" & "+str(int(round(EventCountZero['Wprime_2300_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2300_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2300_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2300_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2300_MixRLWprime'])))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2500 GeV &"+str(EventCountPre['Wprime_2500_MixRLWprime'])+" & "+str(int(round(EventCountZero['Wprime_2500_MixRLWprime'])))+" & "+str(int(round(EventCountOne['Wprime_2500_MixRLWprime'])))+" & "+str(int(round(EventCountGEOne['Wprime_2500_MixRLWprime'])))+" & "+str(int(round(EventCountGETwo['Wprime_2500_MixRLWprime'])))+" & "+str(int(round(EventCountFinal['Wprime_2500_MixRLWprime'])))+"   \\\ \n")
    '''
    txtfile.write("\hline\n")
    txtfile.write("\end{tabular}\n")
    txtfile.write("\\normalsize\n")
    txtfile.write("\end{center}\n")
    txtfile.write("\end{table}\n")
    

    masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2600','2700','2800','2900','3000']

    '''
    print channel
    print 'Double_t eventcount_pre[23] = {'
    for m in masses:
        print str(EventCountPre['Wprime'+m+'Right'])+','
    print '};'
    print 'Double_t eventcount_final[23] = {'
    for m in masses:
        print str(EventCountFinal['Wprime'+m+'Right'])+','
    print '};'
    print 'Double_t signal_LumiXsecOver4[23] = {'
    for m in masses:
        print str(lumi_el*xsec['Wprime'+m+'Right']/4.0)+','
    print '};'
    '''

channel = 'electron'
var = 'elec_1_pt_WprimeCalc'; bin =50; low = 0; high = 7000; 
FillTables(channel, var, bin, low, high) 

channel = 'muon' 
var = 'muon_1_pt_WprimeCalc'; bin =50; low = 0; high = 7000;  
FillTables(channel, var, bin, low, high) 


