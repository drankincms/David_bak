import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack, TFormula
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad

TFormula.SetMaxima(5000,5000,5000)

from array import array

from LoadData_LPC import *

global txtfile
txtfile = open('eventtables_scaleGenTopPt_1and2_withtau_19pt5_fix2btags_withL.txt', 'w') 

CutList = ['preselection',
           'zerobtags',
           'onebtags',
           'twobtags',
           'final1',
           'final2'
    ]

#CutList = ['preselection',
#    ]

def FillTables(channel, varName, bin, low, high):

    VariablesPre = {}
    VariablesZero = {}
    VariablesOne = {}
    VariablesTwo = {}
    VariablesFinal1 = {}
    VariablesFinal2 = {}
 
    WbbHist = {}
    WccHist = {}
  
    EventCountPre = {}
    EventCountZero = {}
    EventCountOne = {}
    EventCountTwo = {}
    EventCountFinal1 = {}
    EventCountFinal2 = {}

    backgroundPre = 0
    backgroundZero = 0
    backgroundOne = 0
    backgroundTwo = 0
    backgroundFinal1 = 0
    backgroundFinal2 = 0
   
    doqcd = 'False'

    if (channel == 'electron'):
        if doqcd == 'True':
            List = ['Data_el','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph','QCD_Pt_80_170_EM','QCD_Pt_170_250_EM','QCD_Pt_250_350_EM','QCD_Pt_350_EM']
        else:
            List = ['Data_el','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']
    if (channel == 'muon'):
        List = ['Data_mu','ZJets_M50','WW', 'WJets', 'T_t', 'Tbar_t', 'T_tW', 'Tbar_tW', 'T_s', 'Tbar_s', 'TTbar_Madgraph']

    SignalList = ['Wprime800Right','Wprime900Right','Wprime1000Right','Wprime1100Right','Wprime1200Right','Wprime1300Right','Wprime1400Right','Wprime1500Right','Wprime1600Right','Wprime1700Right','Wprime1800Right','Wprime1900Right','Wprime2000Right','Wprime2100Right','Wprime2200Right','Wprime2300Right','Wprime2400Right','Wprime2500Right','Wprime2600Right','Wprime2700Right','Wprime2800Right','Wprime2900Right','Wprime3000Right','Wprime1800Left','Wprime800Left','Wprime2000Left','Wprime2500Left','Wprime3000Left']
    List.extend(SignalList)


    for cutlabel in CutList:


        if (channel == 'electron'):      
            cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 50 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20 &&  Muon_DeltaR_LjetsTopoCalcNew > 0.3' 
        if (channel == 'muon'):      
            cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 50 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > 0.3'
        

        cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
        cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
        cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light

        if cutlabel == 'zerobtags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==0) '
        if cutlabel == 'onebtags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==1) '
        if cutlabel == 'twobtags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==2)'
        if cutlabel == 'ge1btags': cut = cut + ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)>=1) '
        if cutlabel == 'ge2btags': cut = cut + ' && ( (jet_0_tag_WprimeCalc==1 && (jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) || (jet_1_tag_WprimeCalc==1 && (jet_0_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) ) '
        if cutlabel == 'final1':    cut = cut + ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) == 1 ) && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 85  && Jet1Jet2_Pt_LjetsTopoCalcNew > 140 '
        if cutlabel == 'final2':    cut = cut + ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) == 2 ) && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 85  && Jet1Jet2_Pt_LjetsTopoCalcNew > 140 '
        if cutlabel == 'finalge1':    cut = cut + ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) >= 1 ) && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > 85  && Jet1Jet2_Pt_LjetsTopoCalcNew > 140 '

        
        SFWjmu = 0.83        ## myHF120lep50, gen level top pt
        SFWcmu = 1.0*1.66   ## myHF120lep50, gen level top pt
        SFWbmu = 1.0*1.21   ## myHF120lep50, gen level top pt                       
        #SFWjmu = 1.0   ## noHF
        #SFWcmu = 1.66  ## noHF
        #SFWbmu = 1.21  ## noHF
              
        weight_ttbar_avg = 'weight_TopPt_WprimeCalc'        
                        
        if (channel=='electron'):
            weight = '( ((0.973*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
        if (channel=='muon'):
            weight = 'weight_PU_ABCD_PileUpCalc*weight_MuonEff_WprimeCalc'

        j = 0
        for Type in List:

            histName = varName+Type+cutlabel+channel
            #print histName
          
            WccHist[cutlabel] = TH1D('WccHist'+histName, 'WccHist'+histName, bin, low, high)
            WbbHist[cutlabel] = TH1D('WbbHist'+histName, 'WbbHist'+histName, bin, low, high)

            SF = 1.0    
            if (Type.startswith('Wprime')): SF = (2.*1.116)/3.

            if cutlabel == 'preselection': VariablesPre[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'zerobtags': VariablesZero[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'onebtags': VariablesOne[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'twobtags': VariablesTwo[Type] = TH1D(histName, histName, bin, low, high)  
            if cutlabel == 'final1': VariablesFinal1[Type] = TH1D(histName, histName, bin, low, high)
            if cutlabel == 'final2': VariablesFinal2[Type] = TH1D(histName, histName, bin, low, high)
             
            #print cut         
            if Type.startswith('Data'):
                Trees[Type].Draw(var + " >> " + histName, "(" + cut + ")", 'goff')
            elif Type == 'WJets':       
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + str(SFWjmu) + ")*(" + cut + cutwjj + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WbbHist"+histName, "("+weight+")*(" + str(SFWbmu) + ")*(" + cut + cutwbb + ")", 'goff')
                Trees[Type].Draw(var + " >> " + "WccHist"+histName, "("+weight+")*(" + str(SFWcmu) + ")*(" + cut + cutwcc + ")", 'goff')
                
                #if cutlabel == 'preselection':
                #    VariablesPre[Type].Add(WbbHist)
                #    VariablesPre[Type].Add(WccHist)
                #if cutlabel == 'zerobtags':
                #    VariablesZero[Type].Add(WbbHist)
                #    VariablesZero[Type].Add(WccHist)
                #if cutlabel == 'onebtags':
                #    VariablesOne[Type].Add(WbbHist)
                #    VariablesOne[Type].Add(WccHist)
                #if cutlabel == 'ge1btags':
                #    VariablesGEOne[Type].Add(WbbHist)
                #    VariablesGEOne[Type].Add(WccHist)
                #if cutlabel == 'ge2btags':
                #    VariablesGETwo[Type].Add(WbbHist)
                #    VariablesGETwo[Type].Add(WccHist)
                #if cutlabel == 'final':
                #    VariablesFinal[Type].Add(WbbHist)
                #    VariablesFinal[Type].Add(WccHist)
            elif Type.startswith('TTbar'):
                Trees[Type].Draw(var + " >> " + histName, "("+weight_ttbar_avg+")*("+weight+")*(" + cut + ")", 'goff')
            else:
                Trees[Type].Draw(var + " >> " + histName, "("+weight+")*(" + cut + ")", 'goff')
                #if (Type.startswith("TTbar")):
                #    Trees[Type].SetScanField(0)    
                #    Trees[Type].Scan('run_CommonCalc:lumi_CommonCalc:event_CommonCalc:corr_met_WprimeCalc:jet_0_pt_WprimeCalc:jet_1_pt_WprimeCalc:elec_1_pt_WprimeCalc:elec_1_eta_WprimeCalc:elec_1_RelIso_WprimeCalc:muon_1_pt_WprimeCalc:muon_1_eta_WprimeCalc:muon_1_RelIso_WprimeCalc', cut)


            if (channel == 'electron'):
                lumi = lumi_el
            if (channel == 'muon'):
                lumi = lumi_mu
            
            if (not Type.startswith('Data')):
                if cutlabel == 'preselection':
                    #print "Preselection "+Type+" Before scaling "+str(VariablesPre[Type].Integral())  
                    VariablesPre[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    #print "Preselection "+Type+" After scaling "+str(VariablesPre[Type].Integral())
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'zerobtags':
                    VariablesZero[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'onebtags':
                    VariablesOne[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'twobtags':
                    VariablesTwo[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'final1':
                    VariablesFinal1[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                if cutlabel == 'final2':
                    VariablesFinal2[Type].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                    if Type.startswith('WJets'): 
                        WbbHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )
                        WccHist[cutlabel].Scale ( (SF*lumi*xsec[Type]/Nevents[Type]) )

 
                if cutlabel == 'preselection':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundPre += VariablesPre[Type].Integral()
                    EventCountPre[Type] = VariablesPre[Type].Integral()
                    #print "EventCounPre "+Type+" "+str(EventCountPre[Type])
                    if Type.startswith('WJets'): 
                        EventCountPre['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountPre['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundPre += WbbHist[cutlabel].Integral()
                        backgroundPre += WccHist[cutlabel].Integral()
                if cutlabel == 'zerobtags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundZero += VariablesZero[Type].Integral()
                    EventCountZero[Type] = VariablesZero[Type].Integral()
                    if Type.startswith('WJets'): 
                        EventCountZero['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountZero['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundZero += WbbHist[cutlabel].Integral()
                        backgroundZero += WccHist[cutlabel].Integral()
                if cutlabel == 'onebtags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundOne += VariablesOne[Type].Integral()
                    EventCountOne[Type] = VariablesOne[Type].Integral()
                    if Type.startswith('WJets'): 
                        EventCountOne['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountOne['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundOne += WbbHist[cutlabel].Integral()
                        backgroundOne += WccHist[cutlabel].Integral()
                if cutlabel == 'twobtags':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundTwo += VariablesTwo[Type].Integral()
                    EventCountTwo[Type] = VariablesTwo[Type].Integral()
                    if Type.startswith('WJets'): 
                        EventCountTwo['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountTwo['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundTwo += WbbHist[cutlabel].Integral()
                        backgroundTwo += WccHist[cutlabel].Integral()
                if cutlabel == 'final1':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundFinal1 += VariablesFinal1[Type].Integral()
                    EventCountFinal1[Type] = VariablesFinal1[Type].Integral()
                    if Type.startswith('WJets'): 
                        EventCountFinal1['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountFinal1['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundFinal1 += WbbHist[cutlabel].Integral()
                        backgroundFinal1 += WccHist[cutlabel].Integral()
                if cutlabel == 'final2':
                    if Type != 'T_s' and Type != 'Tbar_s' and (not Type.startswith('Wprime')): backgroundFinal2 += VariablesFinal2[Type].Integral()
                    EventCountFinal2[Type] = VariablesFinal2[Type].Integral()
                    if Type.startswith('WJets'): 
                        EventCountFinal2['Wbb'] = WbbHist[cutlabel].Integral()
                        EventCountFinal2['Wcc'] = WccHist[cutlabel].Integral()
                        backgroundFinal2 += WbbHist[cutlabel].Integral()
                        backgroundFinal2 += WccHist[cutlabel].Integral()
                     
            if Type.startswith('Data'):
                
                if cutlabel == 'preselection':
                    EventCountPre[Type] = VariablesPre[Type].Integral()
                if cutlabel == 'zerobtags':
                    EventCountZero[Type] = VariablesZero[Type].Integral()
                if cutlabel == 'onebtags':
                    EventCountOne[Type] = VariablesOne[Type].Integral()
                if cutlabel == 'twobtags':
                    EventCountTwo[Type] = VariablesTwo[Type].Integral()
                if cutlabel == 'final1':
                    EventCountFinal1[Type] = VariablesFinal1[Type].Integral()
                if cutlabel == 'final2':
                    EventCountFinal2[Type] = VariablesFinal2[Type].Integral()

            j=j+1

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
      
        EventCountTwo['QCD_Pt_80_170_EM'] += EventCountTwo['QCD_Pt_170_250_EM']
        EventCountTwo['QCD_Pt_80_170_EM'] += EventCountTwo['QCD_Pt_250_350_EM']
        EventCountTwo['QCD_Pt_80_170_EM'] += EventCountTwo['QCD_Pt_350_EM']
      
        EventCountFinal1['QCD_Pt_80_170_EM'] += EventCountFinal1['QCD_Pt_170_250_EM']
        EventCountFinal1['QCD_Pt_80_170_EM'] += EventCountFinal1['QCD_Pt_250_350_EM']
        EventCountFinal1['QCD_Pt_80_170_EM'] += EventCountFinal1['QCD_Pt_350_EM']
      
        EventCountFinal2['QCD_Pt_80_170_EM'] += EventCountFinal2['QCD_Pt_170_250_EM']
        EventCountFinal2['QCD_Pt_80_170_EM'] += EventCountFinal2['QCD_Pt_250_350_EM']
        EventCountFinal2['QCD_Pt_80_170_EM'] += EventCountFinal2['QCD_Pt_350_EM']

    if (channel == 'electron'): suffix = 'el'
    if (channel == 'muon'): suffix = 'mu'

    txtfile.write("\\begin{table}[!h!tb]\n")
    txtfile.write("\\begin{center}\n")
    txtfile.write("\small\n")
    txtfile.write("      \caption{\n")
    txtfile.write("	Number of selected data, and background  events in the "+channel+" channel. \n")
    txtfile.write("For the background samples, the expectation is computed\n")
    txtfile.write("corresponding to an integrated luminosity of $"+str(lumi)+" \pbinv$. ``Final selection'' refers to the additional cuts of $p_{T}^{top}>$ 85, $p_{T}^{jet1,jet2}>$140, and 130$<m_{top}<$210. \n")
    txtfile.write("\label{tab:"+suffix+"_cut_flow_bkg}\n") 
    txtfile.write("        }\n")
    txtfile.write("\\begin{tabular}{l|c|c|c|c|c|c}\n")
    txtfile.write("\hline\n")
    txtfile.write("Process    & \multicolumn{6}{c}{Number of Events} \\\ \hline\hline\n")
    txtfile.write(" & \multicolumn{4}{c|}{Preselection} & \multicolumn{2}{c}{Final Selection} \\\ \n")
    txtfile.write(" & $\geq$ 0 b-tags & = 0 b-tags & = 1 b-tags \n")
    txtfile.write(" & = 2 b-tags & = 1 b-tags  & = 2 b-tags \\\ \hline\n")
    txtfile.write("{\\bf Data} & "+str(round(EventCountPre['Data_'+suffix],1))+" & "+str(round(EventCountZero['Data_'+suffix],1))+" & "+str(round(EventCountOne['Data_'+suffix],1))+" & "+str(round(EventCountTwo['Data_'+suffix],1))+" & "+str(round(EventCountFinal1['Data_'+suffix],1))+" & "+str(round(EventCountFinal2['Data_'+suffix],1))+"   \\\ \hline\n") 
    txtfile.write("\multicolumn{6}{l}{\\bf Background:}     \\\ \hline\n")
    txtfile.write("$t\overline t$ &"+str(round(EventCountPre['TTbar_Madgraph'],1))+" & "+str(round(EventCountZero['TTbar_Madgraph'],1))+" & "+str(round(EventCountOne['TTbar_Madgraph'],1))+" & "+str(round(EventCountTwo['TTbar_Madgraph'],1))+" & "+str(round(EventCountFinal1['TTbar_Madgraph'],1))+" & "+str(round(EventCountFinal2['TTbar_Madgraph'],1))+"  \\\ \n")
    #txtfile.write("$t\overline t$ &"+str(round(EventCountPre['TTbar_Powheg'],1))+" & "+str(round(EventCountZero['TTbar_Powheg'],1))+" & "+str(round(EventCountOne['TTbar_Powheg'],1))+" & "+str(round(EventCountTwo['TTbar_Powheg'],1))+" & "+str(round(EventCountFinal1['TTbar_Powheg'],1))+" & "+str(round(EventCountFinal2['TTbar_Powheg'],1))+"  \\\ \n")
    txtfile.write("$tqb$ &"+str(round(EventCountPre['T_t'],1))+" & "+str(round(EventCountZero['T_t'],1))+" & "+str(round(EventCountOne['T_t'],1))+" & "+str(round(EventCountTwo['T_t'],1))+" & "+str(round(EventCountFinal1['T_t'],1))+"  & "+str(round(EventCountFinal2['T_t'],1))+"  \\\ \n")
    txtfile.write("$\overline{t}qb$ &"+str(round(EventCountPre['Tbar_t'],1))+" & "+str(round(EventCountZero['Tbar_t'],1))+" & "+str(round(EventCountOne['Tbar_t'],1))+" & "+str(round(EventCountTwo['Tbar_t'],1))+" & "+str(round(EventCountFinal1['Tbar_t'],1))+"  & "+str(round(EventCountFinal2['Tbar_t'],1))+"  \\\ \n")
    txtfile.write("$tW$  &"+str(round(EventCountPre['T_tW'],1))+" & "+str(round(EventCountZero['T_tW'],1))+" & "+str(round(EventCountOne['T_tW'],1))+" & "+str(round(EventCountTwo['T_tW'],1))+" & "+str(round(EventCountFinal1['T_tW'],1))+"  & "+str(round(EventCountFinal2['T_tW'],1))+"  \\\ \n")
    txtfile.write("$\overline{t}W$ &"+str(round(EventCountPre['Tbar_tW'],1))+" & "+str(round(EventCountZero['Tbar_tW'],1))+" & "+str(round(EventCountOne['Tbar_tW'],1))+" & "+str(round(EventCountTwo['Tbar_tW'],1))+" & "+str(round(EventCountFinal1['Tbar_tW'],1))+"  & "+str(round(EventCountFinal2['Tbar_tW'],1))+"  \\\ \n")
    txtfile.write("$tb$ & "+str(round(EventCountPre['T_s'],1))+" & "+str(round(EventCountZero['T_s'],1))+" & "+str(round(EventCountOne['T_s'],1))+" & "+str(round(EventCountTwo['T_s'],1))+" & "+str(round(EventCountFinal1['T_s'],1))+" & "+str(round(EventCountFinal2['T_s'],1))+"    \\\ \n")
    txtfile.write("$\overline{t}b$ & "+str(round(EventCountPre['Tbar_s'],1))+" & "+str(round(EventCountZero['Tbar_s'],1))+" & "+str(round(EventCountOne['Tbar_s'],1))+" & "+str(round(EventCountTwo['Tbar_s'],1))+" & "+str(round(EventCountFinal1['Tbar_s'],1))+" & "+str(round(EventCountFinal2['Tbar_s'],1))+"    \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+light jets &"+str(round(EventCountPre['WJets'],1))+" & "+str(round(EventCountZero['WJets'],1))+" & "+str(round(EventCountOne['WJets'],1))+" & "+str(round(EventCountTwo['WJets'],1))+" & "+str(round(EventCountFinal1['WJets'],1))+"  & "+str(round(EventCountFinal2['WJets'],1))+"  \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+c jets &"+str(round(EventCountPre['Wcc'],1))+" & "+str(round(EventCountZero['Wcc'],1))+" & "+str(round(EventCountOne['Wcc'],1))+" & "+str(round(EventCountTwo['Wcc'],1))+" & "+str(round(EventCountFinal1['Wcc'],1))+"  & "+str(round(EventCountFinal2['Wcc'],1))+"  \\\ \n")
    txtfile.write("$W(\\rightarrow)\ell\\nu$+b jets &"+str(round(EventCountPre['Wbb'],1))+" & "+str(round(EventCountZero['Wbb'],1))+" & "+str(round(EventCountOne['Wbb'],1))+" & "+str(round(EventCountTwo['Wbb'],1))+" & "+str(round(EventCountFinal1['Wbb'],1))+"  & "+str(round(EventCountFinal2['Wbb'],1))+"  \\\ \n")
    txtfile.write("$Z/\gamma^*(\\rightarrow\ell\ell$)+jets &"+str(round(EventCountPre['ZJets_M50'],1))+" & "+str(round(EventCountZero['ZJets_M50'],1))+" & "+str(round(EventCountOne['ZJets_M50'],1))+" & "+str(round(EventCountTwo['ZJets_M50'],1))+" & "+str(round(EventCountFinal1['ZJets_M50'],1))+"  & "+str(round(EventCountFinal2['ZJets_M50'],1))+" \\\ \n")
    txtfile.write("$WW$ &"+str(round(EventCountPre['WW'],1))+" & "+str(round(EventCountZero['WW'],1))+" & "+str(round(EventCountOne['WW'],1))+" & "+str(round(EventCountTwo['WW'],1))+" & "+str(round(EventCountFinal1['WW'],1))+"  & "+str(round(EventCountFinal2['WW'],1))+" \\\ \n")
    if (channel=='electron' and doqcd == 'True'):
        txtfile.write("$QCD$ &"+str(round(EventCountPre['QCD_Pt_80_170_EM'],1))+" & "+str(round(EventCountZero['QCD_Pt_80_170_EM'],1))+" & "+str(round(EventCountOne['QCD_Pt_80_170_EM'],1))+" & "+str(round(EventCountTwo['QCD_Pt_80_170_EM'],1))+" & "+str(round(EventCountFinal1['QCD_Pt_80_170_EM'],1))+"  & "+str(round(EventCountFinal2['QCD_Pt_80_170_EM'],1))+" \\\ \n")
    txtfile.write("{\\bf Total Background} &"+str(round(backgroundPre,1))+" & "+str(round(backgroundZero,1))+" & "+str(round(backgroundOne,1))+" & "+str(round(backgroundTwo,1))+" & "+str(round(backgroundFinal1,1))+" & "+str(round(backgroundFinal2,1))+" \\\ \hline\n")
    txtfile.write("{\\bf MC / Data} &"+str(float(round(backgroundPre/EventCountPre['Data_'+suffix], 3)))+" & "+str(float(round(backgroundZero/EventCountZero['Data_'+suffix], 3)))+" & "+str(float(round(backgroundOne/EventCountOne['Data_'+suffix], 3)))+" & "+str(float(round(backgroundTwo/EventCountTwo['Data_'+suffix], 3)))+" & "+str(float(round(backgroundFinal1/EventCountFinal1['Data_'+suffix], 3)))+" & "+str(float(round(backgroundFinal2/EventCountFinal2['Data_'+suffix], 3)))+" \\\ \hline\n")
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
    EventCountTwo['T_s'] += EventCountTwo['Tbar_s'] 
    EventCountFinal1['T_s'] += EventCountFinal1['Tbar_s'] 
    EventCountFinal2['T_s'] += EventCountFinal2['Tbar_s'] 
    

    txtfile.write("\\begin{table}[!h!tb]\n")
    txtfile.write("\\begin{center}\n")
    txtfile.write("\small\n")
    txtfile.write("      \caption{\n")
    txtfile.write("	Number of selected signal events in the "+channel+" channel. \n")
    txtfile.write("For the signal samples, the expectation is computed\n")
    txtfile.write("corresponding to an integrated luminosity of $"+str(lumi)+" \pbinv$. ``Final selection'' refers to the additional cuts of $p_{T}^{top}>$ 85, $p_{T}^{jet1,jet2}>$140, and 130$<m_{top}<$210. \n")
    txtfile.write("\label{tab:"+suffix+"_cut_flow_bkg}\n") 
    txtfile.write("        }\n")
    txtfile.write("\\begin{tabular}{l|c|c|c|c|c|c}\n")
    txtfile.write("\hline\n")
    txtfile.write("Process    & \multicolumn{6}{c}{Number of Events} \\\ \hline\hline\n")
    txtfile.write(" & \multicolumn{4}{c|}{Preselection} & \multicolumn{2}{c}{Final Selection} \\\ \n")
    txtfile.write(" & $\geq$ 0 b-tags & = 0 b-tags & = 1 b-tags \n")
    txtfile.write(" & = 2 b-tags & = 1 b-tags  & = 2 b-tags \\\ \hline\n")  
    txtfile.write("\multicolumn{6}{l}{\\bf Signal:}     \\\ \hline \n")
    txtfile.write("$tb$ & "+str(round(EventCountPre['T_s'],3))+" & "+str(round(EventCountZero['T_s'],3))+" & "+str(round(EventCountOne['T_s'],3))+" & "+str(round(EventCountTwo['T_s'],3))+" & "+str(round(EventCountFinal1['T_s'],3))+" & "+str(round(EventCountFinal2['T_s'],3))+"    \\\ \n")
    txtfile.write("{\underline {\\bf ${\PWpr}_R\\rightarrow tb$}} & &  &  & & \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 800 GeV &"+str(round(EventCountPre['Wprime800Right'],3))+" & "+str(round(EventCountZero['Wprime800Right'],3))+" & "+str(round(EventCountOne['Wprime800Right'],3))+" & "+str(round(EventCountTwo['Wprime800Right'],3))+" & "+str(round(EventCountFinal1['Wprime800Right'],3))+" & "+str(round(EventCountFinal2['Wprime800Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 900 GeV &"+str(round(EventCountPre['Wprime900Right'],3))+" & "+str(round(EventCountZero['Wprime900Right'],3))+" & "+str(round(EventCountOne['Wprime900Right'],3))+" & "+str(round(EventCountTwo['Wprime900Right'],3))+" & "+str(round(EventCountFinal1['Wprime900Right'],3))+" & "+str(round(EventCountFinal2['Wprime900Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1000 GeV &"+str(round(EventCountPre['Wprime1000Right'],3))+" & "+str(round(EventCountZero['Wprime1000Right'],3))+" & "+str(round(EventCountOne['Wprime1000Right'],3))+" & "+str(round(EventCountTwo['Wprime1000Right'],3))+" & "+str(round(EventCountFinal1['Wprime1000Right'],3))+" & "+str(round(EventCountFinal2['Wprime1000Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1100 GeV &"+str(round(EventCountPre['Wprime1100Right'],3))+" & "+str(round(EventCountZero['Wprime1100Right'],3))+" & "+str(round(EventCountOne['Wprime1100Right'],3))+" & "+str(round(EventCountTwo['Wprime1100Right'],3))+" & "+str(round(EventCountFinal1['Wprime1100Right'],3))+" & "+str(round(EventCountFinal2['Wprime1100Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1200 GeV &"+str(round(EventCountPre['Wprime1200Right'],3))+" & "+str(round(EventCountZero['Wprime1200Right'],3))+" & "+str(round(EventCountOne['Wprime1200Right'],3))+" & "+str(round(EventCountTwo['Wprime1200Right'],3))+" & "+str(round(EventCountFinal1['Wprime1200Right'],3))+" & "+str(round(EventCountFinal2['Wprime1200Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1300 GeV &"+str(round(EventCountPre['Wprime1300Right'],3))+" & "+str(round(EventCountZero['Wprime1300Right'],3))+" & "+str(round(EventCountOne['Wprime1300Right'],3))+" & "+str(round(EventCountTwo['Wprime1300Right'],3))+" & "+str(round(EventCountFinal1['Wprime1300Right'],3))+" & "+str(round(EventCountFinal2['Wprime1300Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1400 GeV &"+str(round(EventCountPre['Wprime1400Right'],3))+" & "+str(round(EventCountZero['Wprime1400Right'],3))+" & "+str(round(EventCountOne['Wprime1400Right'],3))+" & "+str(round(EventCountTwo['Wprime1400Right'],3))+" & "+str(round(EventCountFinal1['Wprime1400Right'],3))+" & "+str(round(EventCountFinal2['Wprime1400Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1500 GeV &"+str(round(EventCountPre['Wprime1500Right'],3))+" & "+str(round(EventCountZero['Wprime1500Right'],3))+" & "+str(round(EventCountOne['Wprime1500Right'],3))+" & "+str(round(EventCountTwo['Wprime1500Right'],3))+" & "+str(round(EventCountFinal1['Wprime1500Right'],3))+" & "+str(round(EventCountFinal2['Wprime1500Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1600 GeV &"+str(round(EventCountPre['Wprime1600Right'],3))+" & "+str(round(EventCountZero['Wprime1600Right'],3))+" & "+str(round(EventCountOne['Wprime1600Right'],3))+" & "+str(round(EventCountTwo['Wprime1600Right'],3))+" & "+str(round(EventCountFinal1['Wprime1600Right'],3))+" & "+str(round(EventCountFinal2['Wprime1600Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1700 GeV &"+str(round(EventCountPre['Wprime1700Right'],3))+" & "+str(round(EventCountZero['Wprime1700Right'],3))+" & "+str(round(EventCountOne['Wprime1700Right'],3))+" & "+str(round(EventCountTwo['Wprime1700Right'],3))+" & "+str(round(EventCountFinal1['Wprime1700Right'],3))+" & "+str(round(EventCountFinal2['Wprime1700Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1800 GeV &"+str(round(EventCountPre['Wprime1800Right'],3))+" & "+str(round(EventCountZero['Wprime1800Right'],3))+" & "+str(round(EventCountOne['Wprime1800Right'],3))+" & "+str(round(EventCountTwo['Wprime1800Right'],3))+" & "+str(round(EventCountFinal1['Wprime1800Right'],3))+" & "+str(round(EventCountFinal2['Wprime1800Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 1900 GeV &"+str(round(EventCountPre['Wprime1900Right'],3))+" & "+str(round(EventCountZero['Wprime1900Right'],3))+" & "+str(round(EventCountOne['Wprime1900Right'],3))+" & "+str(round(EventCountTwo['Wprime1900Right'],3))+" & "+str(round(EventCountFinal1['Wprime1900Right'],3))+" & "+str(round(EventCountFinal2['Wprime1900Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2000 GeV &"+str(round(EventCountPre['Wprime2000Right'],3))+" & "+str(round(EventCountZero['Wprime2000Right'],3))+" & "+str(round(EventCountOne['Wprime2000Right'],3))+" & "+str(round(EventCountTwo['Wprime2000Right'],3))+" & "+str(round(EventCountFinal1['Wprime2000Right'],3))+" & "+str(round(EventCountFinal2['Wprime2000Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2100 GeV &"+str(round(EventCountPre['Wprime2100Right'],3))+" & "+str(round(EventCountZero['Wprime2100Right'],3))+" & "+str(round(EventCountOne['Wprime2100Right'],3))+" & "+str(round(EventCountTwo['Wprime2100Right'],3))+" & "+str(round(EventCountFinal1['Wprime2100Right'],3))+" & "+str(round(EventCountFinal2['Wprime2100Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2200 GeV &"+str(round(EventCountPre['Wprime2200Right'],3))+" & "+str(round(EventCountZero['Wprime2200Right'],3))+" & "+str(round(EventCountOne['Wprime2200Right'],3))+" & "+str(round(EventCountTwo['Wprime2200Right'],3))+" & "+str(round(EventCountFinal1['Wprime2200Right'],3))+" & "+str(round(EventCountFinal2['Wprime2200Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2300 GeV &"+str(round(EventCountPre['Wprime2300Right'],3))+" & "+str(round(EventCountZero['Wprime2300Right'],3))+" & "+str(round(EventCountOne['Wprime2300Right'],3))+" & "+str(round(EventCountTwo['Wprime2300Right'],3))+" & "+str(round(EventCountFinal1['Wprime2300Right'],3))+" & "+str(round(EventCountFinal2['Wprime2300Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2500 GeV &"+str(round(EventCountPre['Wprime2500Right'],3))+" & "+str(round(EventCountZero['Wprime2500Right'],3))+" & "+str(round(EventCountOne['Wprime2500Right'],3))+" & "+str(round(EventCountTwo['Wprime2500Right'],3))+" & "+str(round(EventCountFinal1['Wprime2500Right'],3))+" & "+str(round(EventCountFinal2['Wprime2500Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2400 GeV &"+str(round(EventCountPre['Wprime2400Right'],3))+" & "+str(round(EventCountZero['Wprime2400Right'],3))+" & "+str(round(EventCountOne['Wprime2400Right'],3))+" & "+str(round(EventCountTwo['Wprime2400Right'],3))+" & "+str(round(EventCountFinal1['Wprime2400Right'],3))+" & "+str(round(EventCountFinal2['Wprime2400Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2600 GeV &"+str(round(EventCountPre['Wprime2600Right'],3))+" & "+str(round(EventCountZero['Wprime2600Right'],3))+" & "+str(round(EventCountOne['Wprime2600Right'],3))+" & "+str(round(EventCountTwo['Wprime2600Right'],3))+" & "+str(round(EventCountFinal1['Wprime2600Right'],3))+" & "+str(round(EventCountFinal2['Wprime2600Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2700 GeV &"+str(round(EventCountPre['Wprime2700Right'],3))+" & "+str(round(EventCountZero['Wprime2700Right'],3))+" & "+str(round(EventCountOne['Wprime2700Right'],3))+" & "+str(round(EventCountTwo['Wprime2700Right'],3))+" & "+str(round(EventCountFinal1['Wprime2700Right'],3))+" & "+str(round(EventCountFinal2['Wprime2700Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2800 GeV &"+str(round(EventCountPre['Wprime2800Right'],3))+" & "+str(round(EventCountZero['Wprime2800Right'],3))+" & "+str(round(EventCountOne['Wprime2800Right'],3))+" & "+str(round(EventCountTwo['Wprime2800Right'],3))+" & "+str(round(EventCountFinal1['Wprime2800Right'],3))+" & "+str(round(EventCountFinal2['Wprime2800Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 2900 GeV &"+str(round(EventCountPre['Wprime2900Right'],3))+" & "+str(round(EventCountZero['Wprime2900Right'],3))+" & "+str(round(EventCountOne['Wprime2900Right'],3))+" & "+str(round(EventCountTwo['Wprime2900Right'],3))+" & "+str(round(EventCountFinal1['Wprime2900Right'],3))+" & "+str(round(EventCountFinal2['Wprime2900Right'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_R=$) 3000 GeV &"+str(round(EventCountPre['Wprime3000Right'],3))+" & "+str(round(EventCountZero['Wprime3000Right'],3))+" & "+str(round(EventCountOne['Wprime3000Right'],3))+" & "+str(round(EventCountTwo['Wprime3000Right'],3))+" & "+str(round(EventCountFinal1['Wprime3000Right'],3))+" & "+str(round(EventCountFinal2['Wprime3000Right'],3))+"   \\\ \n")

    
    txtfile.write("{\underline {\\bf ${\PWpr}_L\\rightarrow tb$}}  &  & &  & & \\\ \n") 
    txtfile.write("M(${\PWpr}_L=$)800 GeV &"+str(round(EventCountPre['Wprime800Left'],3))+" & "+str(round(EventCountZero['Wprime800Left'],3))+" & "+str(round(EventCountOne['Wprime800Left'],3))+" & "+str(round(EventCountTwo['Wprime800Left'],3))+" & "+str(round(EventCountFinal1['Wprime800Left'],3))+" & "+str(round(EventCountFinal2['Wprime800Left'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_L=$)1800 GeV &"+str(round(EventCountPre['Wprime1800Left'],3))+" & "+str(round(EventCountZero['Wprime1800Left'],3))+" & "+str(round(EventCountOne['Wprime1800Left'],3))+" & "+str(round(EventCountTwo['Wprime1800Left'],3))+" & "+str(round(EventCountFinal1['Wprime1800Left'],3))+" & "+str(round(EventCountFinal2['Wprime1800Left'],3))+"   \\\ \n") 
    txtfile.write("M(${\PWpr}_L=$)2000 GeV &"+str(round(EventCountPre['Wprime2000Left'],3))+" & "+str(round(EventCountZero['Wprime2000Left'],3))+" & "+str(round(EventCountOne['Wprime2000Left'],3))+" & "+str(round(EventCountTwo['Wprime2000Left'],3))+" & "+str(round(EventCountFinal1['Wprime2000Left'],3))+" & "+str(round(EventCountFinal2['Wprime2000Left'],3))+"   \\\ \n") 
    txtfile.write("M(${\PWpr}_L=$)2500 GeV &"+str(round(EventCountPre['Wprime2500Left'],3))+" & "+str(round(EventCountZero['Wprime2500Left'],3))+" & "+str(round(EventCountOne['Wprime2500Left'],3))+" & "+str(round(EventCountTwo['Wprime2500Left'],3))+" & "+str(round(EventCountFinal1['Wprime2500Left'],3))+" & "+str(round(EventCountFinal2['Wprime2500Left'],3))+"   \\\ \n") 
    txtfile.write("M(${\PWpr}_L=$)3000 GeV &"+str(round(EventCountPre['Wprime3000Left'],3))+" & "+str(round(EventCountZero['Wprime3000Left'],3))+" & "+str(round(EventCountOne['Wprime3000Left'],3))+" & "+str(round(EventCountTwo['Wprime3000Left'],3))+" & "+str(round(EventCountFinal1['Wprime3000Left'],3))+" & "+str(round(EventCountFinal2['Wprime3000Left'],3))+"   \\\ \n")
    

    '''
    txtfile.write("{\underline {\\bf ${\PWpr}_{LR}\\rightarrow tb$}}&  & & & & \\\ \n")   
    txtfile.write("M(${\PWpr}_{LR}=$)800 GeV &"+str(round(EventCountPre['Wprime800_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime800_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime800_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime800_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime800_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime800_MixRLWprime'],3))+"    \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)900 GeV &"+str(round(EventCountPre['Wprime900_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime900_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime900_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime900_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime900_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime900_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1000 GeV &"+str(round(EventCountPre['Wprime1000_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1000_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1000_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1000_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1000_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1000_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1100 GeV &"+str(round(EventCountPre['Wprime1100_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1100_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1100_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1100_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1100_MixRLWprime'],3))+"  & "+str(round(EventCountFinal2['Wprime1100_MixRLWprime'],3))+"  \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1200 GeV & "+str(round(EventCountPre['Wprime1200_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1200_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1200_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1200_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1200_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1200_MixRLWprime'],3))+"  \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1300 GeV &"+str(round(EventCountPre['Wprime1300_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1300_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1300_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1300_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1300_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1300_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1400 GeV &"+str(round(EventCountPre['Wprime1400_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1400_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1400_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1400_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1400_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1400_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1500 GeV &"+str(round(EventCountPre['Wprime1500_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1500_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1500_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1500_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1500_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1500_MixRLWprime'],3))+"    \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1600 GeV &"+str(round(EventCountPre['Wprime1600_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1600_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1600_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1600_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1600_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1600_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1700 GeV &"+str(round(EventCountPre['Wprime1700_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1700_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1700_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1700_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1700_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1700_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)1900 GeV &"+str(round(EventCountPre['Wprime1900_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime1900_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime1900_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime1900_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime1900_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime1900_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2100 GeV &"+str(round(EventCountPre['Wprime2100_MixRLWprime'],3))+" & "+str(round(EventCountZero['Wprime2100_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime2100_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime2100_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime2100_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime2100_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2300 GeV &"+str(EventCountPre['Wprime2300_MixRLWprime'])+" & "+str(round(EventCountZero['Wprime2300_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime2300_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime2300_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime2300_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime2300_MixRLWprime'],3))+"   \\\ \n")
    txtfile.write("M(${\PWpr}_{LR}=$)2500 GeV &"+str(EventCountPre['Wprime2500_MixRLWprime'])+" & "+str(round(EventCountZero['Wprime2500_MixRLWprime'],3))+" & "+str(round(EventCountOne['Wprime2500_MixRLWprime'],3))+" & "+str(round(EventCountTwo['Wprime2500_MixRLWprime'],3))+" & "+str(round(EventCountFinal1['Wprime2500_MixRLWprime'],3))+" & "+str(round(EventCountFinal2['Wprime2500_MixRLWprime'],3))+"   \\\ \n")
    '''
    txtfile.write("\hline\n")
    txtfile.write("\end{tabular}\n")
    txtfile.write("\\normalsize\n")
    txtfile.write("\end{center}\n")
    txtfile.write("\end{table}\n")
    

    #masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2600','2700','2800','2900','3000']

    #print channel
    #print 'Double_t eventcount_pre[23] = {'
    #for m in masses:
    #    print str(EventCountPre['Wprime'+m+'Right'])+','
    #print '};'
    #print 'Double_t eventcount_final[23] = {'
    #for m in masses:
    #    print str(EventCountFinal['Wprime'+m+'Right'])+','
    #print '};'
    #print 'Double_t signal_LumiXsecOver4[23] = {'
    #for m in masses:
    #    print str(lumi_el*xsec['Wprime'+m+'Right']/4.0)+','
    #print '};'


channel = 'electron'
var = 'elec_1_pt_WprimeCalc'; bin =50; low = 0; high = 7000; 
FillTables(channel, var, bin, low, high)

channel = 'muon'
var = 'muon_1_pt_WprimeCalc'; bin =50; low = 0; high = 7000; 
FillTables(channel, var, bin, low, high)



