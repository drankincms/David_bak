import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad

from array import array
import copy
import math
import fileinput

#from LoadData import *
from LoadData_LPC import *

List_DataEl = ['Data_el']
List_DataMu = ['Data_mu']

List_Bg = ['WJets','WW','TTbar_Madgraph','ZJets_M50','T_t','Tbar_t','T_tW','Tbar_tW','T_s','Tbar_s',
'WJets_JESUP','WW_JESUP','TTbar_Madgraph_JESUP','ZJets_M50_JESUP','T_t_JESUP','Tbar_t_JESUP','T_tW_JESUP','Tbar_tW_JESUP','T_s_JESUP','Tbar_s_JESUP',
'WJets_JESDOWN','WW_JESDOWN','TTbar_Madgraph_JESDOWN','ZJets_M50_JESDOWN','T_t_JESDOWN','Tbar_t_JESDOWN','T_tW_JESDOWN','Tbar_tW_JESDOWN','T_s_JESDOWN','Tbar_s_JESDOWN',
'WJets_JERUP','WW_JERUP','TTbar_Madgraph_JERUP','ZJets_M50_JERUP','T_t_JERUP','Tbar_t_JERUP','T_tW_JERUP','Tbar_tW_JERUP','T_s_JERUP','Tbar_s_JERUP',
'WJets_JERDOWN','WW_JERDOWN','TTbar_Madgraph_JERDOWN','ZJets_M50_JERDOWN','T_t_JERDOWN','Tbar_t_JERDOWN','T_tW_JERDOWN','Tbar_tW_JERDOWN','T_s_JERDOWN','Tbar_s_JERDOWN',
'WJets_BTAGUP','WW_BTAGUP','TTbar_Madgraph_BTAGUP','ZJets_M50_BTAGUP','T_t_BTAGUP','Tbar_t_BTAGUP','T_tW_BTAGUP','Tbar_tW_BTAGUP','T_s_BTAGUP','Tbar_s_BTAGUP',
'WJets_BTAGDOWN','WW_BTAGDOWN','TTbar_Madgraph_BTAGDOWN','ZJets_M50_BTAGDOWN','T_t_BTAGDOWN','Tbar_t_BTAGDOWN','T_tW_BTAGDOWN','Tbar_tW_BTAGDOWN','T_s_BTAGDOWN','Tbar_s_BTAGDOWN',
'TTbar_Madgraph_MATCHINGUP','TTbar_Madgraph_MATCHINGDOWN','TTbar_Madgraph_SCALEUP','TTbar_Madgraph_SCALEDOWN',
]


List_Right = ['Wprime800Right','Wprime900Right','Wprime1000Right','Wprime1100Right','Wprime1200Right','Wprime1300Right','Wprime1400Right','Wprime1500Right','Wprime1600Right','Wprime1700Right','Wprime1800Right','Wprime1900Right','Wprime2000Right','Wprime2100Right','Wprime2200Right','Wprime2300Right','Wprime2400Right','Wprime2500Right','Wprime2600Right','Wprime2700Right','Wprime2800Right','Wprime2900Right','Wprime3000Right',
'Wprime800Right_JESUP','Wprime900Right_JESUP','Wprime1000Right_JESUP','Wprime1100Right_JESUP','Wprime1200Right_JESUP','Wprime1300Right_JESUP','Wprime1400Right_JESUP','Wprime1500Right_JESUP','Wprime1600Right_JESUP','Wprime1700Right_JESUP','Wprime1800Right_JESUP','Wprime1900Right_JESUP','Wprime2000Right_JESUP','Wprime2100Right_JESUP','Wprime2200Right_JESUP','Wprime2300Right_JESUP','Wprime2400Right_JESUP','Wprime2500Right_JESUP','Wprime2600Right_JESUP','Wprime2700Right_JESUP','Wprime2800Right_JESUP','Wprime2900Right_JESUP','Wprime3000Right_JESUP', 
'Wprime800Right_JESDOWN','Wprime900Right_JESDOWN','Wprime1000Right_JESDOWN','Wprime1100Right_JESDOWN','Wprime1200Right_JESDOWN','Wprime1300Right_JESDOWN','Wprime1400Right_JESDOWN','Wprime1500Right_JESDOWN','Wprime1600Right_JESDOWN','Wprime1700Right_JESDOWN','Wprime1800Right_JESDOWN','Wprime1900Right_JESDOWN','Wprime2000Right_JESDOWN','Wprime2100Right_JESDOWN','Wprime2200Right_JESDOWN','Wprime2300Right_JESDOWN','Wprime2400Right_JESDOWN','Wprime2500Right_JESDOWN','Wprime2600Right_JESDOWN','Wprime2700Right_JESDOWN','Wprime2800Right_JESDOWN','Wprime2900Right_JESDOWN','Wprime3000Right_JESDOWN', 
'Wprime800Right_JERUP','Wprime900Right_JERUP','Wprime1000Right_JERUP','Wprime1100Right_JERUP','Wprime1200Right_JERUP','Wprime1300Right_JERUP','Wprime1400Right_JERUP','Wprime1500Right_JERUP','Wprime1600Right_JERUP','Wprime1700Right_JERUP','Wprime1800Right_JERUP','Wprime1900Right_JERUP','Wprime2000Right_JERUP','Wprime2100Right_JERUP','Wprime2200Right_JERUP','Wprime2300Right_JERUP','Wprime2400Right_JERUP','Wprime2500Right_JERUP','Wprime2600Right_JERUP','Wprime2700Right_JERUP','Wprime2800Right_JERUP','Wprime2900Right_JERUP','Wprime3000Right_JERUP', 
'Wprime800Right_JERDOWN','Wprime900Right_JERDOWN','Wprime1000Right_JERDOWN','Wprime1100Right_JERDOWN','Wprime1200Right_JERDOWN','Wprime1300Right_JERDOWN','Wprime1400Right_JERDOWN','Wprime1500Right_JERDOWN','Wprime1600Right_JERDOWN','Wprime1700Right_JERDOWN','Wprime1800Right_JERDOWN','Wprime1900Right_JERDOWN','Wprime2000Right_JERDOWN','Wprime2100Right_JERDOWN','Wprime2200Right_JERDOWN','Wprime2300Right_JERDOWN','Wprime2400Right_JERDOWN','Wprime2500Right_JERDOWN','Wprime2600Right_JERDOWN','Wprime2700Right_JERDOWN','Wprime2800Right_JERDOWN','Wprime2900Right_JERDOWN','Wprime3000Right_JERDOWN', 
'Wprime800Right_BTAGUP','Wprime900Right_BTAGUP','Wprime1000Right_BTAGUP','Wprime1100Right_BTAGUP','Wprime1200Right_BTAGUP','Wprime1300Right_BTAGUP','Wprime1400Right_BTAGUP','Wprime1500Right_BTAGUP','Wprime1600Right_BTAGUP','Wprime1700Right_BTAGUP','Wprime1800Right_BTAGUP','Wprime1900Right_BTAGUP','Wprime2000Right_BTAGUP','Wprime2100Right_BTAGUP','Wprime2200Right_BTAGUP','Wprime2300Right_BTAGUP','Wprime2400Right_BTAGUP','Wprime2500Right_BTAGUP','Wprime2600Right_BTAGUP','Wprime2700Right_BTAGUP','Wprime2800Right_BTAGUP','Wprime2900Right_BTAGUP','Wprime3000Right_BTAGUP', 
'Wprime800Right_BTAGDOWN','Wprime900Right_BTAGDOWN','Wprime1000Right_BTAGDOWN','Wprime1100Right_BTAGDOWN','Wprime1200Right_BTAGDOWN','Wprime1300Right_BTAGDOWN','Wprime1400Right_BTAGDOWN','Wprime1500Right_BTAGDOWN','Wprime1600Right_BTAGDOWN','Wprime1700Right_BTAGDOWN','Wprime1800Right_BTAGDOWN','Wprime1900Right_BTAGDOWN','Wprime2000Right_BTAGDOWN','Wprime2100Right_BTAGDOWN','Wprime2200Right_BTAGDOWN','Wprime2300Right_BTAGDOWN','Wprime2400Right_BTAGDOWN','Wprime2500Right_BTAGDOWN','Wprime2600Right_BTAGDOWN','Wprime2700Right_BTAGDOWN','Wprime2800Right_BTAGDOWN','Wprime2900Right_BTAGDOWN','Wprime3000Right_BTAGDOWN', 
] 

List_Left = [
'Wprime800Left','Wprime900Left','Wprime1000Left','Wprime1100Left','Wprime1200Left','Wprime1300Left','Wprime1400Left','Wprime1500Left','Wprime1600Left','Wprime1700Left','Wprime1800Left','Wprime1900Left','Wprime2000Left','Wprime2100Left','Wprime2200Left','Wprime2300Left','Wprime2400Left','Wprime2500Left','Wprime2700Left','Wprime2800Left','Wprime2900Left',
'Wprime800Left_JESUP','Wprime900Left_JESUP','Wprime1000Left_JESUP','Wprime1100Left_JESUP','Wprime1200Left_JESUP','Wprime1300Left_JESUP','Wprime1400Left_JESUP','Wprime1500Left_JESUP','Wprime1600Left_JESUP','Wprime1700Left_JESUP','Wprime1800Left_JESUP','Wprime1900Left_JESUP','Wprime2000Left_JESUP','Wprime2100Left_JESUP','Wprime2200Left_JESUP','Wprime2300Left_JESUP','Wprime2400Left_JESUP','Wprime2500Left_JESUP','Wprime2700Left_JESUP','Wprime2800Left_JESUP','Wprime2900Left_JESUP',
'Wprime800Left_JESDOWN','Wprime900Left_JESDOWN','Wprime1000Left_JESDOWN','Wprime1100Left_JESDOWN','Wprime1200Left_JESDOWN','Wprime1300Left_JESDOWN','Wprime1400Left_JESDOWN','Wprime1500Left_JESDOWN','Wprime1600Left_JESDOWN','Wprime1700Left_JESDOWN','Wprime1800Left_JESDOWN','Wprime1900Left_JESDOWN','Wprime2000Left_JESDOWN','Wprime2100Left_JESDOWN','Wprime2200Left_JESDOWN','Wprime2300Left_JESDOWN','Wprime2400Left_JESDOWN','Wprime2500Left_JESDOWN','Wprime2700Left_JESDOWN','Wprime2800Left_JESDOWN','Wprime2900Left_JESDOWN',
'Wprime800Left_JERUP','Wprime900Left_JERUP','Wprime1000Left_JERUP','Wprime1100Left_JERUP','Wprime1200Left_JERUP','Wprime1300Left_JERUP','Wprime1400Left_JERUP','Wprime1500Left_JERUP','Wprime1600Left_JERUP','Wprime1700Left_JERUP','Wprime1800Left_JERUP','Wprime1900Left_JERUP','Wprime2000Left_JERUP','Wprime2100Left_JERUP','Wprime2200Left_JERUP','Wprime2300Left_JERUP','Wprime2400Left_JERUP','Wprime2500Left_JERUP','Wprime2700Left_JERUP','Wprime2800Left_JERUP','Wprime2900Left_JERUP',
'Wprime800Left_JERDOWN','Wprime900Left_JERDOWN','Wprime1000Left_JERDOWN','Wprime1100Left_JERDOWN','Wprime1200Left_JERDOWN','Wprime1300Left_JERDOWN','Wprime1400Left_JERDOWN','Wprime1500Left_JERDOWN','Wprime1600Left_JERDOWN','Wprime1700Left_JERDOWN','Wprime1800Left_JERDOWN','Wprime1900Left_JERDOWN','Wprime2000Left_JERDOWN','Wprime2100Left_JERDOWN','Wprime2200Left_JERDOWN','Wprime2300Left_JERDOWN','Wprime2400Left_JERDOWN','Wprime2500Left_JERDOWN','Wprime2700Left_JERDOWN','Wprime2800Left_JERDOWN','Wprime2900Left_JERDOWN',
'Wprime800Left_BTAGUP','Wprime900Left_BTAGUP','Wprime1000Left_BTAGUP','Wprime1100Left_BTAGUP','Wprime1200Left_BTAGUP','Wprime1300Left_BTAGUP','Wprime1400Left_BTAGUP','Wprime1500Left_BTAGUP','Wprime1600Left_BTAGUP','Wprime1700Left_BTAGUP','Wprime1800Left_BTAGUP','Wprime1900Left_BTAGUP','Wprime2000Left_BTAGUP','Wprime2100Left_BTAGUP','Wprime2200Left_BTAGUP','Wprime2300Left_BTAGUP','Wprime2400Left_BTAGUP','Wprime2500Left_BTAGUP','Wprime2700Left_BTAGUP','Wprime2800Left_BTAGUP','Wprime2900Left_BTAGUP',
'Wprime800Left_BTAGDOWN','Wprime900Left_BTAGDOWN','Wprime1000Left_BTAGDOWN','Wprime1100Left_BTAGDOWN','Wprime1200Left_BTAGDOWN','Wprime1300Left_BTAGDOWN','Wprime1400Left_BTAGDOWN','Wprime1500Left_BTAGDOWN','Wprime1600Left_BTAGDOWN','Wprime1700Left_BTAGDOWN','Wprime1800Left_BTAGDOWN','Wprime1900Left_BTAGDOWN','Wprime2000Left_BTAGDOWN','Wprime2100Left_BTAGDOWN','Wprime2200Left_BTAGDOWN','Wprime2300Left_BTAGDOWN','Wprime2400Left_BTAGDOWN','Wprime2500Left_BTAGDOWN','Wprime2700Left_BTAGDOWN','Wprime2800Left_BTAGDOWN','Wprime2900Left_BTAGDOWN',
]

List_Mix = [
'Wprime800Mix','Wprime900Mix','Wprime1000Mix','Wprime1100Mix','Wprime1200Mix','Wprime1300Mix','Wprime1400Mix','Wprime1500Mix','Wprime1600Mix','Wprime1700Mix','Wprime1800Mix','Wprime1900Mix','Wprime2000Mix','Wprime2100Mix','Wprime2200Mix','Wprime2300Mix','Wprime2400Mix','Wprime2500Mix','Wprime2700Mix','Wprime2800Mix','Wprime2900Mix',
'Wprime800Mix_JESUP','Wprime900Mix_JESUP','Wprime1000Mix_JESUP','Wprime1100Mix_JESUP','Wprime1200Mix_JESUP','Wprime1300Mix_JESUP','Wprime1400Mix_JESUP','Wprime1500Mix_JESUP','Wprime1600Mix_JESUP','Wprime1700Mix_JESUP','Wprime1800Mix_JESUP','Wprime1900Mix_JESUP','Wprime2000Mix_JESUP','Wprime2100Mix_JESUP','Wprime2200Mix_JESUP','Wprime2300Mix_JESUP','Wprime2400Mix_JESUP','Wprime2500Mix_JESUP','Wprime2700Mix_JESUP','Wprime2800Mix_JESUP','Wprime2900Mix_JESUP',
'Wprime800Mix_JESDOWN','Wprime900Mix_JESDOWN','Wprime1000Mix_JESDOWN','Wprime1100Mix_JESDOWN','Wprime1200Mix_JESDOWN','Wprime1300Mix_JESDOWN','Wprime1400Mix_JESDOWN','Wprime1500Mix_JESDOWN','Wprime1600Mix_JESDOWN','Wprime1700Mix_JESDOWN','Wprime1800Mix_JESDOWN','Wprime1900Mix_JESDOWN','Wprime2000Mix_JESDOWN','Wprime2100Mix_JESDOWN','Wprime2200Mix_JESDOWN','Wprime2300Mix_JESDOWN','Wprime2400Mix_JESDOWN','Wprime2500Mix_JESDOWN','Wprime2700Mix_JESDOWN','Wprime2800Mix_JESDOWN','Wprime2900Mix_JESDOWN',
'Wprime800Mix_JERUP','Wprime900Mix_JERUP','Wprime1000Mix_JERUP','Wprime1100Mix_JERUP','Wprime1200Mix_JERUP','Wprime1300Mix_JERUP','Wprime1400Mix_JERUP','Wprime1500Mix_JERUP','Wprime1600Mix_JERUP','Wprime1700Mix_JERUP','Wprime1800Mix_JERUP','Wprime1900Mix_JERUP','Wprime2000Mix_JERUP','Wprime2100Mix_JERUP','Wprime2200Mix_JERUP','Wprime2300Mix_JERUP','Wprime2400Mix_JERUP','Wprime2500Mix_JERUP','Wprime2700Mix_JERUP','Wprime2800Mix_JERUP','Wprime2900Mix_JERUP',
'Wprime800Mix_JERDOWN','Wprime900Mix_JERDOWN','Wprime1000Mix_JERDOWN','Wprime1100Mix_JERDOWN','Wprime1200Mix_JERDOWN','Wprime1300Mix_JERDOWN','Wprime1400Mix_JERDOWN','Wprime1500Mix_JERDOWN','Wprime1600Mix_JERDOWN','Wprime1700Mix_JERDOWN','Wprime1800Mix_JERDOWN','Wprime1900Mix_JERDOWN','Wprime2000Mix_JERDOWN','Wprime2100Mix_JERDOWN','Wprime2200Mix_JERDOWN','Wprime2300Mix_JERDOWN','Wprime2400Mix_JERDOWN','Wprime2500Mix_JERDOWN','Wprime2700Mix_JERDOWN','Wprime2800Mix_JERDOWN','Wprime2900Mix_JERDOWN',
'Wprime800Mix_BTAGUP','Wprime900Mix_BTAGUP','Wprime1000Mix_BTAGUP','Wprime1100Mix_BTAGUP','Wprime1200Mix_BTAGUP','Wprime1300Mix_BTAGUP','Wprime1400Mix_BTAGUP','Wprime1500Mix_BTAGUP','Wprime1600Mix_BTAGUP','Wprime1700Mix_BTAGUP','Wprime1800Mix_BTAGUP','Wprime1900Mix_BTAGUP','Wprime2000Mix_BTAGUP','Wprime2100Mix_BTAGUP','Wprime2200Mix_BTAGUP','Wprime2300Mix_BTAGUP','Wprime2400Mix_BTAGUP','Wprime2500Mix_BTAGUP','Wprime2700Mix_BTAGUP','Wprime2800Mix_BTAGUP','Wprime2900Mix_BTAGUP',
'Wprime800Mix_BTAGDOWN','Wprime900Mix_BTAGDOWN','Wprime1000Mix_BTAGDOWN','Wprime1100Mix_BTAGDOWN','Wprime1200Mix_BTAGDOWN','Wprime1300Mix_BTAGDOWN','Wprime1400Mix_BTAGDOWN','Wprime1500Mix_BTAGDOWN','Wprime1600Mix_BTAGDOWN','Wprime1700Mix_BTAGDOWN','Wprime1800Mix_BTAGDOWN','Wprime1900Mix_BTAGDOWN','Wprime2000Mix_BTAGDOWN','Wprime2100Mix_BTAGDOWN','Wprime2200Mix_BTAGDOWN','Wprime2300Mix_BTAGDOWN','Wprime2400Mix_BTAGDOWN','Wprime2500Mix_BTAGDOWN','Wprime2700Mix_BTAGDOWN','Wprime2800Mix_BTAGDOWN','Wprime2900Mix_BTAGDOWN',
]


def create_yields_list(channel,varName, toppt, j1j2pt, lepjetdR, bin, low, high, ylabel, xlabel, save, wprime, btags, List_to_use):

    if (channel == 'electron'): ch = '_el'
    if (channel == 'muon'): ch = '_mu'

    doQCD = True
    if (channel == 'electron' and doQCD):
        List_to_use.extend(['QCD_Pt_80_170_EM','QCD_Pt_170_250_EM','QCD_Pt_250_350_EM','QCD_Pt_350_EM',]]])

    pyfile = open("yields_For2DLimits_06Jun_finalbins_ScaleGenTopPt_QCD/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_dr03_lep50.py","w") 
        

    doTTbarWeight = 'True'

    if (channel == 'electron'):      
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 50 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > '+lepjetdR 

    if (channel == 'muon'):      
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 50 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > '+lepjetdR 
   
    print varName
                 

    if btags == 'zerobtags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==0) '
    if btags == 'onebtags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==1) '
    if btags == 'ge1btags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)>=1) '
    if btags == 'twobtags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==2) '
    if btags == 'ge2btags': cutbtag = ' && ( (jet_0_tag_WprimeCalc==1 && (jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) || (jet_1_tag_WprimeCalc==1 && (jet_0_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) ) '

    if btags == 'final': cutbtag = ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) >= 1 )'
    if btags == 'final': cut = cut + ' && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > '+toppt+' && Jet1Jet2_Pt_LjetsTopoCalcNew > '+j1j2pt


    cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
    cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
    cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light


    if (doTTbarWeight == 'True'):
        SFWjmu = 0.82        ## myHF120lep50, scale ttbar
        SFWcmu = 1.0*1.66   ## myHF120lep50, scale ttbar
        SFWbmu = 1.0*1.21   ## myHF120lep50, scale ttbar

        SFWjmuPlus = 0.82*0.87        ## myHF120lep50, scale ttbar
        SFWcmuPlus = 1.0*1.66*1.15   ## myHF120lep50, scale ttbar
        SFWbmuPlus = 1.0*1.21*1.15   ## myHF120lep50, scale ttbar
            
        SFWjmuMinus = 0.82*1.13        ## myHF120lep50, scale ttbar
        SFWcmuMinus = 1.0*1.66*0.85   ## myHF120lep50, scale ttbar
        SFWbmuMinus = 1.0*1.21*0.85   ## myHF120lep50, scale ttbar

        #weight_ttbarminus = '(1.188-0.001124*(BestTop_Pt_LjetsTopoCalcNew))'
        #weight_ttbarplus = '1.0'
        #weight_ttbar = '(5.809*TMath::Landau(BestTop_Pt_LjetsTopoCalcNew,149.1,154))'   
        weight_ttbarminus = 'weight_TopPt_WprimeCalc'
        weight_ttbarplus = '1.0'
        weight_ttbar = 'weight_TopPt_WprimeCalc'
                        
    if (doTTbarWeight == 'False'):
        SFWjmu = 0.82        ## myHF120lep50
        SFWcmu = 0.93*1.66   ## myHF120lep50
        SFWbmu = 0.93*1.21   ## myHF120lep50

        SFWjmuPlus = 0.82*0.87        ## myHF120lep50
        SFWcmuPlus = 0.93*1.66*1.15   ## myHF120lep50
        SFWbmuPlus = 0.93*1.21*1.15   ## myHF120lep50
            
        SFWjmuMinus = 0.82*1.13        ## myHF120lep50
        SFWcmuMinus = 0.93*1.66*0.85   ## myHF120lep50
        SFWbmuMinus = 0.93*1.21*0.85   ## myHF120lep50
        
        weight_ttbarminus = '1.0'
        weight_ttbarplus = '1.0'
        weight_ttbar = '1.0'   
    
    cutzerobtags = ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc+jet_2_tag_WprimeCalc+jet_3_tag_WprimeCalc+jet_4_tag_WprimeCalc+jet_5_tag_WprimeCalc+jet_6_tag_WprimeCalc+jet_7_tag_WprimeCalc+jet_8_tag_WprimeCalc+jet_9_tag_WprimeCalc) == 0 )' 
 
    print wprime 
    #print cut + cutbtag
 
    List = List_to_use
    
    Variables = {}
    VariablesPUup = {}
    Variables0tag = {}
    VariablesHFup = {}
    VariablesHFdown = {}
    VariablesTTbarShapeUp = {}
    VariablesTTbarShapeDown = {}  

    background = 0
    j = 0

    nominalwprime = 'False'

    for Type in List:
    
        if (channel == 'electron'):
            prefix = 'elec_invmass_' + btags + '__'
        if (channel == 'muon'):
            prefix = 'mu_invmass_' + btags + '__'

        suffix = ''
        
        if (Type=='Data_el' or Type=='Data_mu'): suffix = 'DATA' + Type
        if (Type=='WJets'): suffix = 'wjets' + Type
        if (Type=='WW'): suffix = 'scaledntb' + Type
        if (Type=='ZJets_M50' or Type=='WW' or Type=='T_t' or Type=='Tbar_t' or Type=='T_tW' or Type=='Tbar_tW'): suffix = 'scaledntb' + Type
        if (Type=='TTbar_Madgraph'): suffix = 'ttbar' + Type
        if (Type=='T_s' or Type=='Tbar_s'): suffix = 'tb' + Type

        if (Type== 'WJets_JESUP'): suffix = 'wjets_jesUp' + Type
        if (Type=='TTbar_Madgraph_JESUP' or Type=='WW_JESUP' or Type=='ZJets_M50_JESUP' or Type=='T_t_JESUP' or Type=='Tbar_t_JESUP' or Type=='T_tW_JESUP' or Type=='Tbar_tW_JESUP'): suffix = 'scaledntb_jesUp' + Type
        if (Type=='TTbar_Madgraph_JESUP'): suffix = 'scaledall_jesUp' + Type
        if (Type=='T_s_JESUP' or Type=='Tbar_s_JESUP'): suffix = 'tb_jesUp' + Type

        if (Type== 'WJets_JESDOWN'): suffix = 'wjets_jesDown' + Type
        if (Type=='TTbar_Madgraph_JESDOWN' or Type=='WW_JESDOWN' or Type=='ZJets_M50_JESDOWN' or Type=='T_t_JESDOWN' or Type=='Tbar_t_JESDOWN' or Type=='T_tW_JESDOWN' or Type=='Tbar_tW_JESDOWN'): suffix = 'scaledntb_jesDown' + Type
        if (Type=='TTbar_Madgraph_JESDOWN'): suffix = 'scaledall_jesDown' + Type
        if (Type=='T_s_JESDOWN' or Type=='Tbar_s_JESDOWN'): suffix = 'tb_jesDown' + Type

        if (Type== 'WJets_JERUP'): suffix = 'wjets_jerUp' + Type
        if (Type=='TTbar_Madgraph_JERUP' or Type=='WW_JERUP' or Type=='ZJets_M50_JERUP' or Type=='T_t_JERUP' or Type=='Tbar_t_JERUP' or Type=='T_tW_JERUP' or Type=='Tbar_tW_JERUP'): suffix = 'scaledntb_jerUp' + Type
        if (Type=='TTbar_Madgraph_JERUP'): suffix = 'scaledall_jerUp' + Type
        if (Type=='T_s_JERUP' or Type=='Tbar_s_JERUP'): suffix = 'tb_jerUp' + Type

        if (Type== 'WJets_JERDOWN'): suffix = 'wjets_jerDown' + Type
        if (Type=='TTbar_Madgraph_JERDOWN' or Type=='WW_JERDOWN' or Type=='ZJets_M50_JERDOWN' or Type=='T_t_JERDOWN' or Type=='Tbar_t_JERDOWN' or Type=='T_tW_JERDOWN' or Type=='Tbar_tW_JERDOWN'): suffix = 'scaledntb_jerDown' + Type
        if (Type=='TTbar_Madgraph_JERDOWN'): suffix = 'scaledall_jerDown' + Type
        if (Type=='T_s_JERDOWN' or Type=='Tbar_s_JERDOWN'): suffix = 'tb_jerDown' + Type

        if (Type== 'WJets_BTAGUP'): suffix = 'wjets_btagUp' + Type
        if (Type=='TTbar_Madgraph_BTAGUP' or Type=='WW_BTAGUP' or Type=='ZJets_M50_BTAGUP' or Type=='T_t_BTAGUP' or Type=='Tbar_t_BTAGUP' or Type=='T_tW_BTAGUP' or Type=='Tbar_tW_BTAGUP'): suffix = 'scaledntb_btagUp' + Type
        if (Type=='TTbar_Madgraph_BTAGUP'): suffix = 'scaledall_btagUp' + Type
        if (Type=='T_s_BTAGUP' or Type=='Tbar_s_BTAGUP'): suffix = 'tb_btagUp' + Type

        if (Type== 'WJets_BTAGDOWN'): suffix = 'wjets_btagDown' + Type
        if (Type=='TTbar_Madgraph_BTAGDOWN' or Type=='WW_BTAGDOWN' or Type=='ZJets_M50_BTAGDOWN' or Type=='T_t_BTAGDOWN' or Type=='Tbar_t_BTAGDOWN' or Type=='T_tW_BTAGDOWN' or Type=='Tbar_tW_BTAGDOWN'): suffix = 'scaledntb_btagDown' + Type
        if (Type=='TTbar_Madgraph_BTAGDOWN'): suffix = 'scaledall_btagDown' + Type
        if (Type=='T_s_BTAGDOWN' or Type=='Tbar_s_BTAGDOWN'): suffix = 'tb_btagDown' + Type


        if (Type == 'WJets_SCALEUP'): suffix = 'wjets_q2scaleUp'
        if (Type == 'WJets_SCALEDOWN'): suffix = 'wjets_q2scaleDown'
        if (Type == 'WJets_MATCHINGUP'): suffix = 'wjets_matchingUp'
        if (Type == 'WJets_MATCHINGDOWN'): suffix = 'wjets_matchingDown'
        if (Type == 'TTbar_Madgraph_SCALEUP'): suffix = 'ttbar__q2scale__plus'
        if (Type == 'TTbar_Madgraph_SCALEDOWN'): suffix = 'ttbar__q2scale__mins'
        if (Type == 'TTbar_Madgraph_MATCHINGUP'): suffix = 'ttbar__matching__plus'
        if (Type == 'TTbar_Madgraph_MATCHINGDOWN'): suffix = 'ttbar__matching__minus'

        w_suffix = 'wp'

        if (wprime != 'ModRight'):
            if (Type.startswith('Wprime800' + wprime)): suffix = w_suffix+'800'
            if (Type.startswith('Wprime900' + wprime)): suffix = w_suffix+'900'
            if (Type.startswith('Wprime1000' + wprime)): suffix = w_suffix+'1000'
            if (Type.startswith('Wprime1100' + wprime)): suffix = w_suffix+'1100'
            if (Type.startswith('Wprime1200' + wprime)): suffix = w_suffix+'1200'
            if (Type.startswith('Wprime1300' + wprime)): suffix = w_suffix+'1300'
            if (Type.startswith('Wprime1400' + wprime)): suffix = w_suffix+'1400'
            if (Type.startswith('Wprime1500' + wprime)): suffix = w_suffix+'1500'
            if (Type.startswith('Wprime1600' + wprime)): suffix = w_suffix+'1600'
            if (Type.startswith('Wprime1700' + wprime)): suffix = w_suffix+'1700'
            if (Type.startswith('Wprime1800' + wprime)): suffix = w_suffix+'1800'
            if (Type.startswith('Wprime1900' + wprime)): suffix = w_suffix+'1900'
            if (Type.startswith('Wprime2000' + wprime)): suffix = w_suffix+'2000'
            if (Type.startswith('Wprime2100' + wprime)): suffix = w_suffix+'2100'
            if (Type.startswith('Wprime2200' + wprime)): suffix = w_suffix+'2200'
            if (Type.startswith('Wprime2300' + wprime)): suffix = w_suffix+'2300'
            if (Type.startswith('Wprime2400' + wprime)): suffix = w_suffix+'2400'
            if (Type.startswith('Wprime2500' + wprime)): suffix = w_suffix+'2500'
            if (Type.startswith('Wprime2600' + wprime)): suffix = w_suffix+'2600'
            if (Type.startswith('Wprime2700' + wprime)): suffix = w_suffix+'2700'
            if (Type.startswith('Wprime2800' + wprime)): suffix = w_suffix+'2800'
            if (Type.startswith('Wprime2900' + wprime)): suffix = w_suffix+'2900'
            if (Type.startswith('Wprime3000' + wprime)): suffix = w_suffix+'3000'

            if (Type.endswith('_JESUP')): suffix += '__jes__plus'
            if (Type.endswith('_JESDOWN')): suffix += '__jes__minus'
            if (Type.endswith('_JERUP')): suffix += '__jer__plus'
            if (Type.endswith('_JERDOWN')): suffix += '__jer__minus'
            if (Type.endswith('_BTAGUP')): suffix += '__btag__plus'
            if (Type.endswith('_BTAGDOWN')): suffix += '__btag__minus'

        else:

            if (Type.startswith('Wprime800Right')): suffix = w_suffix+'800'
            if (Type.startswith('Wprime900Right')): suffix = w_suffix+'900'
            if (Type.startswith('Wprime1000Right')): suffix = w_suffix+'1000'
            if (Type.startswith('Wprime1100Right')): suffix = w_suffix+'1100'
            if (Type.startswith('Wprime1200Right')): suffix = w_suffix+'1200'
            if (Type.startswith('Wprime1300Right')): suffix = w_suffix+'1300'
            if (Type.startswith('Wprime1400Right')): suffix = w_suffix+'1400'
            if (Type.startswith('Wprime1500Right')): suffix = w_suffix+'1500'
            if (Type.startswith('Wprime1600Right')): suffix = w_suffix+'1600'
            if (Type.startswith('Wprime1700Right')): suffix = w_suffix+'1700'
            if (Type.startswith('Wprime1800Right')): suffix = w_suffix+'1800'
            if (Type.startswith('Wprime1900Right')): suffix = w_suffix+'1900'
            if (Type.startswith('Wprime2000Right')): suffix = w_suffix+'2000'
            if (Type.startswith('Wprime2100Right')): suffix = w_suffix+'2100'
            if (Type.startswith('Wprime2200Right')): suffix = w_suffix+'2200'
            if (Type.startswith('Wprime2300Right')): suffix = w_suffix+'2300'
            if (Type.startswith('Wprime2400Right')): suffix = w_suffix+'2400'
            if (Type.startswith('Wprime2500Right')): suffix = w_suffix+'2500'
            if (Type.startswith('Wprime2600Right')): suffix = w_suffix+'2600'
            if (Type.startswith('Wprime2700Right')): suffix = w_suffix+'2700'
            if (Type.startswith('Wprime2800Right')): suffix = w_suffix+'2800'
            if (Type.startswith('Wprime2900Right')): suffix = w_suffix+'2900'
            if (Type.startswith('Wprime3000Right')): suffix = w_suffix+'3000'

            if (Type.endswith('_JESUP')): suffix += '__jes__plus'
            if (Type.endswith('_JESDOWN')): suffix += '__jes__minus'
            if (Type.endswith('_JERUP')): suffix += '__jer__plus'
            if (Type.endswith('_JERDOWN')): suffix += '__jer__minus'
            if (Type.endswith('_BTAGUP')): suffix += '__btag__plus'
            if (Type.endswith('_BTAGDOWN')): suffix += '__btag__minus'

      
        histName = prefix+suffix+'varbin'
        histNamePUup = prefix+suffix+'varbin'+'__PU__plus'
        histNamePUdown = prefix+suffix+'varbin'+'__PU__minus'
        histName0tag = prefix+suffix+'varbin'+'__0tag__plus'
        histNameHFup = prefix+suffix+'varbin'+'__hf__plus'
        histNameHFdown = prefix+suffix+'varbin'+'__hf__minus'
        histNameTTbarShapeUp = prefix+suffix+'varbin'+'__ttbarshape__plus'
        histNameTTbarShapeDown = prefix+suffix+'varbin'+'__ttbarshape__minus'

        histNamePre = prefix+suffix+'varbin'+'Pre'


        Variables[Type] = TH1D(histName, histName, bin, array('d',xlow))  
        Variables[Type].Sumw2()

        if (channel=='electron'):
            weight = '( ((0.973*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
            weightPUup = '( ((0.973*weight_PU_ABCD735_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD735_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
            SF = 1.0
        if (channel=='muon'):
            weight = 'weight_PU_ABCD_PileUpCalc*weight_MuonEff_WprimeCalc'
            weightPUup = 'weight_PU_ABCD735_PileUpCalc*weight_MuonEff_WprimeCalc'
            SF = 1.0
                
        if (Type == 'WJets'):
            WccHist = TH1D('WccHist', 'WccHist', bin,array('d',xlow))
            WbbHist = TH1D('WbbHist', 'WbbHist', bin,array('d',xlow))
            WccHist.Sumw2()
            WbbHist.Sumw2()
            WccHistPUup = TH1D('WccHistPUup', 'WccHistPUup', bin,array('d',xlow))
            WbbHistPUup = TH1D('WbbHistPUup', 'WbbHistPUup', bin,array('d',xlow))
            WccHistPUup.Sumw2()
            WbbHistPUup.Sumw2()
            WccHistHFup = TH1D('WccHistHFup', 'WccHistHFup', bin,array('d',xlow))
            WbbHistHFup = TH1D('WbbHistHFup', 'WbbHistHFup', bin,array('d',xlow))
            WccHistHFup.Sumw2()
            WbbHistHFup.Sumw2()
            WccHistHFdown = TH1D('WccHistHFdown', 'WccHistHFdown', bin,array('d',xlow))
            WbbHistHFdown = TH1D('WbbHistHFdown', 'WbbHistHFdown', bin,array('d',xlow))
            WccHistHFdown.Sumw2()
            WbbHistHFdown.Sumw2()
            WccHistPre = TH1D('WccHistPre', 'WccHistPre', bin,array('d',xlow))
            WbbHistPre = TH1D('WbbHistPre', 'WbbHistPre', bin,array('d',xlow))
            WccHistPre.Sumw2()
            WbbHistPre.Sumw2()

        if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):                      
            VariablesPUup[Type] = TH1D(histNamePUup, histNamePUup, bin, array('d',xlow))
            #VariablesPUdown[Type] = TH1D(histNamePUdown, histNamePUdown, bin, array('d',xlow))
            Variables0tag[Type] = TH1D(histName0tag, histName0tag, bin, array('d',xlow))
            VariablesHFup[Type] = TH1D(histNameHFup, histNameHFup, bin, array('d',xlow))
            VariablesHFdown[Type] = TH1D(histNameHFdown, histNameHFdown, bin, array('d',xlow))
            VariablesTTbarShapeUp[Type] = TH1D(histNameTTbarShapeUp, histNameTTbarShapeUp, bin, array('d',xlow))
            VariablesTTbarShapeDown[Type] = TH1D(histNameTTbarShapeDown, histNameTTbarShapeDown, bin, array('d',xlow))

            VariablesPUup[Type].Sumw2()
            #VariablesPUdown[Type].Sumw2()
            Variables0tag[Type].Sumw2()
            VariablesHFup[Type].Sumw2()
            VariablesHFdown[Type].Sumw2()
            VariablesTTbarShapeUp[Type].Sumw2()
            VariablesTTbarShapeDown[Type].Sumw2()


        #print Type
        if (Type.startswith('Data')):
            Trees[Type].Draw(var + " >> " + histName, "(" + cut + cutbtag + ")", 'goff')
            # 0 tag for data-driven shape
            Trees[Type].Draw(var + " >> " + histName0tag, "(" + cut + cutzerobtags + ")", 'goff')
        elif (Type.startswith('WJets')): 

            Trees[Type].Draw(var+" >> "+histName,"("+weight+")*("+str(SFWjmu)+")*("+cut+cutwjj+cutbtag+")",'goff')
            Trees[Type].Draw(var+" >> "+"WbbHist","("+weight+")*("+str(SFWbmu)+")*("+cut+cutwbb+cutbtag+")",'goff')
            Trees[Type].Draw(var+" >> "+"WccHist","("+weight+")*("+str(SFWcmu)+")*("+cut+cutwcc+cutbtag+")",'goff') 
            Variables[Type].Add(WbbHist)
            Variables[Type].Add(WccHist)
            if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ): 
                # Pile Up 
                Trees[Type].Draw(var+" >> "+histNamePUup,"("+weightPUup+")*("+str(SFWjmu)+")*("+cut+cutwjj+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WbbHistPUup","("+weightPUup+")*("+str(SFWbmu)+")*("+cut+cutwbb+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WccHistPUup","("+weightPUup+")*("+str(SFWcmu)+")*("+cut+cutwcc+cutbtag+")",'goff') 
                VariablesPUup[Type].Add(WbbHistPUup)
                VariablesPUup[Type].Add(WccHistPUup)
                # 0 tag shape (one sided, so only one histogram)
                Trees[Type].Draw(var+" >> "+histName0tag,"("+weight+")*("+cut+cutzerobtags+")",'goff')
                # H.F. k-factor  
                Trees[Type].Draw(var+" >> "+histNameHFup,"("+weight+")*("+str(SFWjmuPlus)+")*("+cut+cutwjj+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WbbHistHFup","("+weight+")*("+str(SFWbmuPlus)+")*("+cut+cutwbb+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WccHistHFup","("+weight+")*("+str(SFWcmuPlus)+")*("+cut+cutwcc+cutbtag+")",'goff')
                VariablesHFup[Type].Add(WbbHistHFup)
                VariablesHFup[Type].Add(WccHistHFup)
                Trees[Type].Draw(var+" >> "+histNameHFdown,"("+weight+")*("+str(SFWjmuMinus)+")*("+cut+cutwjj+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WbbHistHFdown","("+weight+")*("+str(SFWbmuMinus)+")*("+cut+cutwbb+cutbtag+")",'goff')
                Trees[Type].Draw(var+" >> "+"WccHistHFdown","("+weight+")*("+str(SFWcmuMinus)+")*("+cut+cutwcc+cutbtag+")",'goff')
                VariablesHFdown[Type].Add(WbbHistHFdown)
                VariablesHFdown[Type].Add(WccHistHFdown)
        elif (not Type.startswith('T')):
            Trees[Type].Draw(var + " >> " + histNamePre, "("+weight+")*("+cut+")", 'goff')
            Trees[Type].Draw(var + " >> " + histName, "("+weight+")*("+cut+cutbtag+")", 'goff')
            if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):                      
                Trees[Type].Draw(var + " >> " + histNamePUup, "("+weightPUup+")*("+cut+cutbtag+")", 'goff')
                #Trees[Type].Draw(var + " >> " + histNamePUdown, "("+weightPUdown+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histName0tag, "("+weight+")*("+cut+cutzerobtags+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFup, "("+weight+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFdown, "("+weight+")*("+cut+cutbtag+")", 'goff')
        elif (Type.startswith('TTbar')):
            Trees[Type].Draw(var + " >> " + histName, "("+weight+")*("+weight_ttbar+")*("+cut+cutbtag+")", 'goff')
            if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):
                Trees[Type].Draw(var + " >> " + histNamePUup, "("+weightPUup+")*("+weight_ttbar+")*("+cut+cutbtag+")", 'goff')
                #Trees[Type].Draw(var + " >> " + histNamePUdown, "("+weightPUdown+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histName0tag, "("+weight+")*("+weight_ttbar+")*("+cut+cutzerobtags+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFup, "("+weight+")*("+weight_ttbar+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFdown, "("+weight+")*("+weight_ttbar+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameTTbarShapeUp, "("+weight+")*("+weight_ttbarplus+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameTTbarShapeDown, "("+weight+")*("+weight_ttbarminus+")*("+cut+cutbtag+")", 'goff')
        else:
            #print 'cut = ',cut+cutbtag
            Trees[Type].Draw(var + " >> " + histName, "("+weight+")*("+cut+cutbtag+")", 'goff')
            if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):                      
                Trees[Type].Draw(var + " >> " + histNamePUup, "("+weightPUup+")*("+cut+cutbtag+")", 'goff')
                #Trees[Type].Draw(var + " >> " + histNamePUdown, "("+weightPUdown+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histName0tag, "("+weight+")*("+cut+cutzerobtags+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFup, "("+weight+")*("+cut+cutbtag+")", 'goff')
                Trees[Type].Draw(var + " >> " + histNameHFdown, "("+weight+")*("+cut+cutbtag+")", 'goff')

        if (not Type.startswith('Data')):
            #print 'EVENTS Before Scaling FOR ',Type,' = ',Variables[Type].Integral()
            #print str(SF),' ',str(lumi),' ',str(xsec[Type]),' ',str(Nevents[Type])
                                    
            if (channel=='electron'): lumi = lumi_el
            if (channel=='muon'): lumi = lumi_mu

            if Variables[Type].Integral() != 0:
                Variables[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )

                pyfile.write("    yields[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(Variables[Type].Integral())+"\n")  
                     
                if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):
                    
                    VariablesPUup[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    #VariablesPUdown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )
                    Variables0tag[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )  
                    VariablesHFup[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    VariablesHFdown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )
                    VariablesTTbarShapeUp[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    VariablesTTbarShapeDown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )

                    pyfile.write("    yieldsPUup[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(VariablesPUup[Type].Integral())+"\n") 
                    pyfile.write("    yields0tag[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(Variables0tag[Type].Integral())+"\n")  
                    pyfile.write("    yieldsHFup[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(VariablesHFup[Type].Integral())+"\n") 
                    pyfile.write("    yieldsHFdown[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(VariablesHFdown[Type].Integral())+"\n") 
                    pyfile.write("    yieldsTTbarShapeUp[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(VariablesTTbarShapeUp[Type].Integral())+"\n") 
                    pyfile.write("    yieldsTTbarShapeDown[mass+coup+'_"+Type+ch+"_"+btags+"'] = "+str(VariablesTTbarShapeDown[Type].Integral())+"\n") 
                                

    pyfile.close()
  
    for line in fileinput.input(pyfile.name, inplace=1):
        if '_el_el' in line:
            line = line.replace('_el_el','_el')
        if '_mu_mu' in line:
            line = line.replace('_mu_mu','_mu')
        if 'Data' in line:
            line = line.replace('Data','data')
        if 'WJets' in line:
            line = line.replace('WJets','wjets')
        if 'ZJets_M50' in line:
            line = line.replace('ZJets_M50','zjets')
        if 'WW' in line:
            line = line.replace('WW','ww')
        if 'T_t' in line:
            line = line.replace('T_t','t')
        if 'Tbar_t' in line:
            line = line.replace('Tbar_t','bt')
        if 'tW' in line:
            line = line.replace('tW','tw')
        if 'T_s' in line:
            line = line.replace('T_s','s')
        if 'Tbar_s' in line:
            line = line.replace('Tbar_s','bs')
        if 'TTbar_Madgraph' in line:
            line = line.replace('TTbar_Madgraph','ttbar')
        #if 'QCD_Pt_80_170_EM' in line:
        #    line = line.replace('QCD_Pt_80_170_EM','qcd80to170EM')
        #if 'QCD_Pt_170_250_EM' in line:
        #    line = line.replace('QCD_Pt_170_250_EM','qcd170to250EM')
        #if 'QCD_Pt_250_350_EM' in line:
        #    line = line.replace('QCD_Pt_250_350_EM','qcd250to350EM')
        #if 'QCD_Pt_350_EM' in line:
        #    line = line.replace('QCD_Pt_350_EM','qcd350EM')
        if 'Wprime' in line:
            line = line.replace('Wprime','wp')
        if 'Right' in line:
            line = line.replace('Right','R')
        if 'Left' in line:
            line = line.replace('Left','L')
        if 'Mix' in line:
            line = line.replace('Mix','RL')
        print line    

  
#wprime = 'Right'
var = 'BestJetJet2W_M_LjetsTopoCalcNew'; high = 4000; xaxis = "W' invariant mass [GeV/c^{2}]"; yaxis = 'Events / 10 GeV'; save = 'BestJetJet2W_M'

#btags = 0
#btags = 'final'
#btags = 'ge1btags'

xlow  = [ 100, 300,  400, 500, 600, 700, 800, 900,  1000, 1100, 1200,  1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 4000]
bins = len(xlow)-1 

List_DataBgEl_RightWprime = copy.copy(List_DataEl) 
List_DataBgEl_RightWprime.extend(List_Bg) 
List_DataBgEl_RightWprime.extend(List_Right) 

List_DataBgMu_RightWprime = copy.copy(List_DataMu) 
List_DataBgMu_RightWprime.extend(List_Bg) 
List_DataBgMu_RightWprime.extend(List_Right) 


wprime = 'Right'
btags = 'ge1btags'
channel = 'electron'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_RightWprime)
channel = 'muon'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_RightWprime)
btags = 'onebtags'
channel = 'electron'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_RightWprime) 
channel = 'muon'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_RightWprime) 
btags = 'twobtags'
channel = 'electron'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_RightWprime)
channel = 'muon'
create_yields_list(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_RightWprime)

