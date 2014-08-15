import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad

from array import array
import copy
import math

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

List_ModRight = ['Wprime800Right','Wprime900Right','Wprime1000Right','Wprime1100Right','Wprime1200Right','Wprime1300Right','Wprime1400Right','Wprime1500Right','Wprime1600Right','Wprime1700Right','Wprime1800Right','Wprime1900Right','Wprime2000Right','Wprime2100Right','Wprime2200Right','Wprime2300Right','Wprime2400Right','Wprime2500Right','Wprime2700Right','Wprime2800Right','Wprime2900Right',
'Wprime800Right_JESUP','Wprime900Right_JESUP','Wprime1000Right_JESUP','Wprime1100Right_JESUP','Wprime1200Right_JESUP','Wprime1300Right_JESUP','Wprime1400Right_JESUP','Wprime1500Right_JESUP','Wprime1600Right_JESUP','Wprime1700Right_JESUP','Wprime1800Right_JESUP','Wprime1900Right_JESUP','Wprime2000Right_JESUP','Wprime2100Right_JESUP','Wprime2200Right_JESUP','Wprime2300Right_JESUP','Wprime2400Right_JESUP','Wprime2500Right_JESUP','Wprime2700Right_JESUP','Wprime2800Right_JESUP','Wprime2900Right_JESUP',
'Wprime800Right_JESDOWN','Wprime900Right_JESDOWN','Wprime1000Right_JESDOWN','Wprime1100Right_JESDOWN','Wprime1200Right_JESDOWN','Wprime1300Right_JESDOWN','Wprime1400Right_JESDOWN','Wprime1500Right_JESDOWN','Wprime1600Right_JESDOWN','Wprime1700Right_JESDOWN','Wprime1800Right_JESDOWN','Wprime1900Right_JESDOWN','Wprime2000Right_JESDOWN','Wprime2100Right_JESDOWN','Wprime2200Right_JESDOWN','Wprime2300Right_JESDOWN','Wprime2400Right_JESDOWN','Wprime2500Right_JESDOWN','Wprime2700Right_JESDOWN','Wprime2800Right_JESDOWN','Wprime2900Right_JESDOWN',
'Wprime800Right_JERUP','Wprime900Right_JERUP','Wprime1000Right_JERUP','Wprime1100Right_JERUP','Wprime1200Right_JERUP','Wprime1300Right_JERUP','Wprime1400Right_JERUP','Wprime1500Right_JERUP','Wprime1600Right_JERUP','Wprime1700Right_JERUP','Wprime1800Right_JERUP','Wprime1900Right_JERUP','Wprime2000Right_JERUP','Wprime2100Right_JERUP','Wprime2200Right_JERUP','Wprime2300Right_JERUP','Wprime2400Right_JERUP','Wprime2500Right_JERUP','Wprime2700Right_JERUP','Wprime2800Right_JERUP','Wprime2900Right_JERUP',
'Wprime800Right_JERDOWN','Wprime900Right_JERDOWN','Wprime1000Right_JERDOWN','Wprime1100Right_JERDOWN','Wprime1200Right_JERDOWN','Wprime1300Right_JERDOWN','Wprime1400Right_JERDOWN','Wprime1500Right_JERDOWN','Wprime1600Right_JERDOWN','Wprime1700Right_JERDOWN','Wprime1800Right_JERDOWN','Wprime1900Right_JERDOWN','Wprime2000Right_JERDOWN','Wprime2100Right_JERDOWN','Wprime2200Right_JERDOWN','Wprime2300Right_JERDOWN','Wprime2400Right_JERDOWN','Wprime2500Right_JERDOWN','Wprime2700Right_JERDOWN','Wprime2800Right_JERDOWN','Wprime2900Right_JERDOWN',
'Wprime800Right_BTAGUP','Wprime900Right_BTAGUP','Wprime1000Right_BTAGUP','Wprime1100Right_BTAGUP','Wprime1200Right_BTAGUP','Wprime1300Right_BTAGUP','Wprime1400Right_BTAGUP','Wprime1500Right_BTAGUP','Wprime1600Right_BTAGUP','Wprime1700Right_BTAGUP','Wprime1800Right_BTAGUP','Wprime1900Right_BTAGUP','Wprime2000Right_BTAGUP','Wprime2100Right_BTAGUP','Wprime2200Right_BTAGUP','Wprime2300Right_BTAGUP','Wprime2400Right_BTAGUP','Wprime2500Right_BTAGUP','Wprime2700Right_BTAGUP','Wprime2800Right_BTAGUP','Wprime2900Right_BTAGUP',
'Wprime800Right_BTAGDOWN','Wprime900Right_BTAGDOWN','Wprime1000Right_BTAGDOWN','Wprime1100Right_BTAGDOWN','Wprime1200Right_BTAGDOWN','Wprime1300Right_BTAGDOWN','Wprime1400Right_BTAGDOWN','Wprime1500Right_BTAGDOWN','Wprime1600Right_BTAGDOWN','Wprime1700Right_BTAGDOWN','Wprime1800Right_BTAGDOWN','Wprime1900Right_BTAGDOWN','Wprime2000Right_BTAGDOWN','Wprime2100Right_BTAGDOWN','Wprime2200Right_BTAGDOWN','Wprime2300Right_BTAGDOWN','Wprime2400Right_BTAGDOWN','Wprime2500Right_BTAGDOWN','Wprime2700Right_BTAGDOWN','Wprime2800Right_BTAGDOWN','Wprime2900Right_BTAGDOWN',
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


def plot_DataVsMc(channel,varName, toppt, j1j2pt, lepjetdR, bin, low, high, ylabel, xlabel, save, wprime, btags, List_to_use):


    doTTbarWeight = 'True'

    if (channel == 'electron'):      
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && elec_1_pt_WprimeCalc > 50 && abs(elec_1_eta_WprimeCalc) < 2.5 && elec_1_RelIso_WprimeCalc < 0.1 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > '+lepjetdR 

    if (channel == 'muon'):      
        cut = 'jet_0_pt_WprimeCalc >= 120 && jet_1_pt_WprimeCalc >= 40 && muon_1_pt_WprimeCalc > 50 && abs(muon_1_eta_WprimeCalc) < 2.1 && muon_1_RelIso_WprimeCalc < 0.12 && corr_met_WprimeCalc > 20 && Muon_DeltaR_LjetsTopoCalcNew > '+lepjetdR 
   
    print varName
                 

    if btags == 'zerobtags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==0) '
    if btags == 'onebtags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)==1) '
    if btags == 'ge1btags': cutbtag = ' && ((jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc)>=1) '
    if btags == 'ge2btags': cutbtag = ' && ( (jet_0_tag_WprimeCalc==1 && (jet_1_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) || (jet_1_tag_WprimeCalc==1 && (jet_0_tag_WprimeCalc + jet_2_tag_WprimeCalc + jet_3_tag_WprimeCalc + jet_4_tag_WprimeCalc + jet_5_tag_WprimeCalc + jet_6_tag_WprimeCalc + jet_7_tag_WprimeCalc + jet_8_tag_WprimeCalc + jet_9_tag_WprimeCalc) >= 1 ) ) '
    if btags == 'final': cutbtag = ' && ( (jet_0_tag_WprimeCalc+jet_1_tag_WprimeCalc ) >= 1 )'
    
    if btags == 'final': cut = cut + ' && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > '+toppt+' && Jet1Jet2_Pt_LjetsTopoCalcNew > '+j1j2pt

    #if btags == 'final': cut = cut + ' && BestTop_Pt_LjetsTopoCalcNew > '+toppt+' && Jet1Jet2_Pt_LjetsTopoCalcNew > '+j1j2pt
    #if btags == 'final': cut = cut + ' && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 && Jet1Jet2_Pt_LjetsTopoCalcNew > '+j1j2pt
    #if btags == 'final': cut = cut + ' && BestTop_LjetsTopoCalcNew > 130 && BestTop_LjetsTopoCalcNew < 210 &&  BestTop_Pt_LjetsTopoCalcNew > '+toppt

    
    cutwbb = ' && n_Bjets_WprimeCalc > 0' # Wb(b)
    cutwcc = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc>0' # Wc(c)
    cutwjj = ' && n_Bjets_WprimeCalc==0 && n_Cjets_WprimeCalc==0' # W+light

    #SFWjmu = 1.08*0.85      ## myHF120
    #SFWcmu = 1.06*0.92*1.66  ## myHF120
    #SFWbmu = 1.06*0.92*1.21  ## myHF120

    #SFWjmuPlus = 1.08*0.85*0.8
    #SFWcmuPlus = 1.06*0.92*1.66*1.27
    #SFWbmuPlus = 1.06*0.92*1.21*1.27

    #SFWjmuMinus = 1.08*0.85*1.2
    #SFWcmuMinus = 1.06*0.92*1.66*0.73
    #SFWbmuMinus = 1.06*0.92*1.21*0.73

    #SFWjmu = 0.86      ## myHF120
    #SFWcmu = 0.95*1.66  ## myHF120
    #SFWbmu = 0.95*1.21  ## myHF120

    #SFWjmuPlus = 0.86*0.8
    #SFWcmuPlus = 0.95*1.66*1.27
    #SFWbmuPlus = 0.95*1.21*1.27

    #SFWjmuMinus = 0.86*1.2
    #SFWcmuMinus = 0.95*1.66*0.73
    #SFWbmuMinus = 0.95*1.21*0.73

    #SFWjmu = 0.82      ## myHF120
    #SFWcmu = 0.93*1.66  ## myHF120
    #SFWbmu = 0.93*1.21  ## myHF120

    #SFWjmuPlus = 0.82*0.87
    #SFWcmuPlus = 0.93*1.66*1.15
    #SFWbmuPlus = 0.93*1.21*1.15

    #SFWjmuMinus = 0.82*1.13
    #SFWcmuMinus = 0.93*1.66*0.85
    #SFWbmuMinus = 0.93*1.21*0.85

    #SFWjmu = 0.93  ## myHF120, no SHyFT
    #SFWcmu = 1.06  ## myHF120, no SHyFT
    #SFWbmu = 1.06  ## myHF120,no SHyFT
            
    #SFWjmuPlus = 0.93*0.8  
    #SFWcmuPlus = 1.06*1.27  
    #SFWbmuPlus = 1.06*1.27  

    #SFWjmuMinus = 0.93*1.2
    #SFWcmuMinus = 1.06*0.73
    #SFWbmuMinus = 1.06*0.73

    SFWjmu = 0.82        ## myHF120lep50, scale ttbar
    SFWcmu = 0.98*1.66   ## myHF120lep50, scale ttbar
    SFWbmu = 0.98*1.21   ## myHF120lep50, scale ttbar

    SFWjmuPlus = 0.82*0.87        ## myHF120lep50, scale ttbar
    SFWcmuPlus = 0.98*1.66*1.15   ## myHF120lep50, scale ttbar
    SFWbmuPlus = 0.98*1.21*1.15   ## myHF120lep50, scale ttbar
            
    SFWjmuMinus = 0.82*1.13        ## myHF120lep50, scale ttbar
    SFWcmuMinus = 0.98*1.66*0.85   ## myHF120lep50, scale ttbar
    SFWbmuMinus = 0.98*1.21*0.85   ## myHF120lep50, scale ttbar


    #weight_ttbarminus = '(0.9822-0.00009837*(BestJetJet2W_M_LjetsTopoCalcNew))'
    #weight_ttbarplus = '(1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew))'
    #weight_ttbar = '(0.5*((1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew))+(0.9822-0.00009837*(BestJetJet2W_M_LjetsTopoCalcNew))))'

    #weight_ttbarminus = '1.0'
    #weight_ttbarplus = '(1.0/(1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew)))'
    #weight_ttbar = '(1.168-0.001209*(BestTop_Pt_LjetsTopoCalcNew))'

    #weight_ttbarminus = '(0.5*((1.098-0.0008829*(BestTop_Pt_LjetsTopoCalcNew))+(1.3-0.001659*(BestTop_Pt_LjetsTopoCalcNew))))'
    #weight_ttbarplus = '1.0'
    #weight_ttbar = '(0.5*((1.098-0.0008829*(BestTop_Pt_LjetsTopoCalcNew))+(1.3-0.001659*(BestTop_Pt_LjetsTopoCalcNew))))'

    weight_ttbarminus = '(1.188-0.001124*(BestTop_Pt_LjetsTopoCalcNew))'
    weight_ttbarplus = '1.0'
    weight_ttbar = '(5.809*TMath::Landau(BestTop_Pt_LjetsTopoCalcNew,149.1,154))'              
            
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
    VariablesPre = {}

    VariablesSmooth = {}
    VariablesSmoothPUup = {}
    VariablesSmooth0tag = {}
    VariablesSmoothHFup = {}
    VariablesSmoothHFdown = {}
    VariablesSmoothTTbarShapeUp = {}
    VariablesSmoothTTbarShapeDown = {} 
    VariablesSmoothPre = {}

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

        histNameSmooth = prefix+suffix
        histNameSmoothPUup = prefix+suffix+'__PU__plus'
        histNameSmoothPUdown = prefix+suffix+'__PU_minus'
        histNameSmooth0tag = prefix+suffix+'0tag'
        histNameSmoothHFup = prefix+suffix+'__hf__plus'
        histNameSmoothHFdown = prefix+suffix+'__hf__minus'
        histNameSmoothTTbarShapeUp = prefix+suffix+'__ttbarshape__plus'
        histNameSmoothTTbarShapeDown = prefix+suffix+'__ttbarshape__minus'

        histNameSmoothPre = prefix+suffix+'Pre'

        Variables[Type] = TH1D(histName, histName, bin, array('d',xlow))  
        Variables[Type].Sumw2()
        VariablesPre[Type] = TH1D(histNamePre, histNamePre, bin, array('d',xlow))  
        VariablesPre[Type].Sumw2()
        VariablesSmooth[Type] = TH1D(histNameSmooth, histNameSmooth, bin, array('d',xlow) ) 
        VariablesSmooth[Type].Sumw2()
        VariablesSmoothPre[Type] = TH1D(histNameSmoothPre, histNameSmoothPre, bin, array('d',xlow) ) 
        VariablesSmoothPre[Type].Sumw2()

        if (channel=='electron'):
            #weight = 'weight_PU_ABCD_PileUpCalc*weight_ElectronEff_WprimeCalc'
            #weightPUup = 'weight_PU_ABCD735_PileUpCalc*weight_ElectronEff_53x_WprimeCalc'
            weight = '( ((0.973*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
            weightPUup = '( ((0.973*weight_PU_ABCD735_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)<1.5)) + ((1.02*weight_PU_ABCD735_PileUpCalc*weight_ElectronEff_53x_WprimeCalc)*(abs(elec_1_eta_WprimeCalc)>1.5 && abs(elec_1_eta_WprimeCalc)<2.5)) )'
            SF = 1.0
        if (channel=='muon'):
            weight = 'weight_PU_ABCD_PileUpCalc*weight_MuonEff_WprimeCalc'
            weightPUup = 'weight_PU_ABCD735_PileUpCalc*weight_MuonEff_WprimeCalc'
            SF = 1.0
        
        #if (Type.startswith('TTbar')): SF = 0.95
        
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
            VariablesSmoothPUup[Type] = TH1D(histNameSmoothPUup, histNameSmoothPUup, bin, array('d',xlow) )
            #VariablesSmoothPUdown[Type] = TH1D(histNameSmoothPUdown, histNameSmoothPUdown, bin, array('d',xlow) )
            VariablesSmooth0tag[Type] = TH1D(histNameSmooth0tag, histNameSmooth0tag, bin, array('d',xlow) )
            VariablesSmoothHFup[Type] = TH1D(histNameSmoothHFup, histNameSmoothHFup, bin, array('d',xlow) )
            VariablesSmoothHFdown[Type] = TH1D(histNameSmoothHFdown, histNameSmoothHFdown, bin, array('d',xlow) )
            VariablesSmoothTTbarShapeUp[Type] = TH1D(histNameSmoothTTbarShapeUp, histNameSmoothTTbarShapeUp, bin, array('d',xlow) )
            VariablesSmoothTTbarShapeDown[Type] = TH1D(histNameSmoothTTbarShapeDown, histNameSmoothTTbarShapeDown, bin, array('d',xlow) )

            VariablesSmoothPUup[Type].Sumw2()
            #VariablesSmoothPUdown[Type].Sumw2()
            VariablesSmooth0tag[Type].Sumw2()
            VariablesSmoothHFup[Type].Sumw2()
            VariablesSmoothHFdown[Type].Sumw2()   
            VariablesSmoothTTbarShapeUp[Type].Sumw2()
            VariablesSmoothTTbarShapeDown[Type].Sumw2()   


        #print Type
        if (Type.startswith('Data')):
            Trees[Type].Draw(var + " >> " + histName, "(" + cut + cutbtag + ")", 'goff')
            # 0 tag for data-driven shape
            Trees[Type].Draw(var + " >> " + histName0tag, "(" + cut + cutzerobtags + ")", 'goff')
        elif (Type.startswith('WJets')): 
            ############################################
            ##### Pretag, to be scaled down later
            ############################################
            # here keep weight_WxsecNoLight_comb in case we need to change the scale factors later             
            Trees[Type].Draw(var+" >> "+histNamePre,"("+weight+")*("+str(SFWjmu)+")*("+cut+cutwjj+")",'goff')
            Trees[Type].Draw(var+" >> "+"WbbHistPre","("+weight+")*("+str(SFWbmu)+")*("+cut+cutwbb+")",'goff')
            Trees[Type].Draw(var+" >> "+"WccHistPre","("+weight+")*("+str(SFWcmu)+")*("+cut+cutwcc+")",'goff')
            #print 'WJets light pre: ',VariablesPre[Type].Integral() 
            VariablesPre[Type].Add(WbbHistPre)
            #print 'WJets light+bb pre: ',VariablesPre[Type].Integral() 
            VariablesPre[Type].Add(WccHistPre)
            #print 'WJets light+bb+cc pre: ',VariablesPre[Type].Integral() 
            ############################################
            ##### >=1, to get the ratio on the fly
            ############################################
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
            #print 'Pre Events before scaling for ',Type,' = ',VariablesPre[Type].Integral()
            #print str(SF),' ',str(lumi),' ',str(xsec_norm[Type]),' ',str(Nevents[Type])
                                    
            if (channel=='electron'): lumi = lumi_el
            if (channel=='muon'): lumi = lumi_mu

            if Variables[Type].Integral() != 0:
                Variables[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )
                
                #if (doTTbarSF == 'True' and Type.startswith('TTbar') ):
                #    f1 = TF1("f1","pol1")
                #    f1.SetParameter(0,0.9964)
                #    f1.SetParameter(1,-0.0001151)
                #    for i in range(Variables[Type].GetNbinsX()):
                #        x = Variables[Type].GetBinCenter(i)
                #        sfi = f1.Eval(x)
                #        Variables[Type].SetBinContent(i,Variables[Type].GetBinContent(i)*sfi)
                #        print Type,' ',x,' old: ',Variables[Type].GetBinContent(i),' new: ',Variables[Type].GetBinContent(i)*sfi
                        
                
                if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):
                    
                    VariablesPUup[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    #VariablesPUdown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )
                    Variables0tag[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )  
                    VariablesHFup[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    VariablesHFdown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )
                    VariablesTTbarShapeUp[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
                    VariablesTTbarShapeDown[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) )

                VariablesPre[Type].Scale ( (SF*lumi*xsec_norm[Type]/Nevents[Type]) ) 
   

        ############################################## 
        ### REBIN TO HAVE EVEN SIZED BINS  
        ### If VariablesSmooth* have the same binning
        ### as Variables*, this does nothing
        ##############################################
       
        for x in range(0,bin+1):       
            VariablesSmooth[Type].SetBinContent(x,Variables[Type].GetBinContent(x) )
            VariablesSmooth[Type].SetBinError(x,Variables[Type].GetBinError(x) )
            VariablesSmoothPre[Type].SetBinContent(x,VariablesPre[Type].GetBinContent(x) )
            VariablesSmoothPre[Type].SetBinError(x,VariablesPre[Type].GetBinError(x) )      
            if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):
                VariablesSmoothPUup[Type].SetBinContent(x,VariablesPUup[Type].GetBinContent(x) )
                VariablesSmoothPUup[Type].SetBinError(x,VariablesPUup[Type].GetBinError(x) )
                #VariablesSmoothPUdown[Type].SetBinContent(x,VariablesPUdown[Type].GetBinContent(x) )
                #VariablesSmoothPUdown[Type].SetBinError(x,VariablesPUdown[Type].GetBinError(x) )
                VariablesSmooth0tag[Type].SetBinContent(x,Variables0tag[Type].GetBinContent(x) )
                VariablesSmoothHFup[Type].SetBinContent(x,VariablesHFup[Type].GetBinContent(x)  )
                VariablesSmoothHFdown[Type].SetBinContent(x,VariablesHFdown[Type].GetBinContent(x) )
                VariablesSmoothHFup[Type].SetBinError(x,VariablesHFup[Type].GetBinError(x)  )
                VariablesSmoothHFdown[Type].SetBinError(x,VariablesHFdown[Type].GetBinError(x) )
                VariablesSmoothTTbarShapeUp[Type].SetBinContent(x,VariablesTTbarShapeUp[Type].GetBinContent(x)  )
                VariablesSmoothTTbarShapeDown[Type].SetBinContent(x,VariablesTTbarShapeDown[Type].GetBinContent(x) )
                VariablesSmoothTTbarShapeUp[Type].SetBinError(x,VariablesTTbarShapeUp[Type].GetBinError(x)  )
                VariablesSmoothTTbarShapeDown[Type].SetBinError(x,VariablesTTbarShapeDown[Type].GetBinError(x) )

        VariablesSmooth[Type].SetEntries(Variables[Type].GetEntries() )
        VariablesSmoothPre[Type].SetEntries(VariablesPre[Type].GetEntries() )
        if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ):
            VariablesSmoothPUup[Type].SetEntries(VariablesPUup[Type].GetEntries() )
            #VariablesSmoothPUdown[Type].SetEntries(VariablesPUdown[Type].GetEntries() )
            VariablesSmooth0tag[Type].SetEntries(Variables0tag[Type].GetEntries() )
            VariablesSmoothHFup[Type].SetEntries(VariablesHFup[Type].GetEntries() )
            VariablesSmoothHFdown[Type].SetEntries(VariablesHFdown[Type].GetEntries() )
            VariablesSmoothTTbarShapeUp[Type].SetEntries(VariablesTTbarShapeUp[Type].GetEntries() )
            VariablesSmoothTTbarShapeDown[Type].SetEntries(VariablesTTbarShapeDown[Type].GetEntries() )

        if (Type.startswith('Data')):
            print 'EVENTS FOR  Data  = ',int(VariablesSmooth[Type].Integral())
            if VariablesSmooth[Type].GetBinContent(bin+1)!=0: print 'OVERFLOW!!!!!!'
            if VariablesSmooth[Type].GetBinContent(0)!=0: print 'UNDERFLOW!!!!!!'

        ######################################### 
        ### Set 0 B.G. Bins to something !=0 
        #########################################

        if (not (Type.startswith('Data'))):
            if ( not Type.startswith('Wprime') ):
                for x in range(1,bin+1):
                    if (VariablesSmooth[Type].GetBinContent(x) < 0.000001 ): 
                        print 'Setting ',VariablesSmooth[Type].GetBinContent(x),' to 10E-6 for bin ',x,' of ',Type
                        VariablesSmooth[Type].SetBinContent(x,0.000001)
                    if ( (not Type.endswith('UP')) and (not Type.endswith('DOWN')) ): 
                        if (VariablesSmoothPUup[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothPUup[Type].SetBinContent(x,0.00001)
                        #if (VariablesSmoothPUdown[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothPUdown[Type].SetBinContent(x,0.00001)
                        if (VariablesSmooth0tag[Type].GetBinContent(x) < 0.000001 ): VariablesSmooth0tag[Type].SetBinContent(x,0.00001)
                        if (VariablesSmoothHFup[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothHFup[Type].SetBinContent(x,0.00001)
                        if (VariablesSmoothHFdown[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothHFdown[Type].SetBinContent(x,0.00001)
                        if (VariablesSmoothTTbarShapeUp[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothTTbarShapeUp[Type].SetBinContent(x,0.00001)
                        if (VariablesSmoothTTbarShapeDown[Type].GetBinContent(x) < 0.000001 ): VariablesSmoothTTbarShapeDown[Type].SetBinContent(x,0.00001)

            print 'SCALED EVENTS FOR ',Type,'  = ',Variables[Type].Integral()
            if ( (not Type.startswith('T')) and (not Type.startswith('Wprime')) ):
                print 'Pre EVENTS FOR ',Type,'  = ',VariablesPre[Type].Integral()          

            if ( (not Type.startswith('Data')) and (not Type.startswith('Wprime')) and (not Type.endswith('UP')) and (not Type.endswith('Down'))):
                background = background + Variables[Type].Integral()


    if (channel == 'electron'):
        chan = 'elec_'
        VariablesSmooth['Data_el'].SetName(chan+"invmass_"+btags+"__DATA")
        VariablesSmooth['Data_el'].SetTitle(chan+"invmass_"+btags+"__DATA")
        VariablesSmooth['Data_el'].Write()
    if (channel == 'muon'):
        chan = 'mu_'
        VariablesSmooth['Data_mu'].SetName(chan+"invmass_"+btags+"__DATA")
        VariablesSmooth['Data_mu'].SetTitle(chan+"invmass_"+btags+"__DATA")
        VariablesSmooth['Data_mu'].Write()

    VariablesSmooth['WW'].Add(VariablesSmooth['ZJets_M50'])
    VariablesSmooth['WJets'].Add(VariablesSmooth['WW'])
    VariablesSmoothPre['WW'].Add(VariablesSmoothPre['ZJets_M50'])
    VariablesSmoothPre['WJets'].Add(VariablesSmoothPre['WW'])
    VariablesSmoothPre['WJets'].SetName(chan+"invmass_"+btags+"__wjets")
    VariablesSmoothPre['WJets'].SetTitle(chan+"invmass_"+btags+"__wjets")

    ######################################### 
    ### RESCALE W+Jets    
    #########################################
    print 'WJets ge1/pre: ',VariablesSmooth['WJets'].Integral()/VariablesSmoothPre['WJets'].Integral()
    VariablesSmoothPre['WJets'].Scale( VariablesSmooth['WJets'].Integral()/VariablesSmoothPre['WJets'].Integral() )
    VariablesSmoothPre['WJets'].Write()
    VariablesSmooth['WJets'].SetName(chan+"invmass_"+btags+"__ge1bwjets")
    VariablesSmooth['WJets'].SetTitle(chan+"invmass_"+btags+"__ge1bwjets")
    #VariablesSmooth['WJets'].Write()

    VariablesSmooth['TTbar_Madgraph'].SetName(chan+"invmass_"+btags+"__ttbar")
    VariablesSmooth['TTbar_Madgraph'].SetTitle(chan+"invmass_"+btags+"__ttbar")
    VariablesSmooth['T_t'].Add(VariablesSmooth['Tbar_t'])
    VariablesSmooth['T_t'].Add(VariablesSmooth['T_tW'])
    VariablesSmooth['T_t'].Add(VariablesSmooth['Tbar_tW'])
    VariablesSmooth['TTbar_Madgraph'].Add( VariablesSmooth['T_t'])

    VariablesSmooth['T_s'].Add(VariablesSmooth['Tbar_s']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph'].Write()
            VariablesSmooth['T_s'].SetName(chan+"invmass_"+btags+"__topstb")
            VariablesSmooth['T_s'].SetTitle(chan+"invmass_"+btags+"__topstb")
            VariablesSmooth['T_s'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph'].Add(VariablesSmooth['T_s'])
            VariablesSmooth['TTbar_Madgraph'].Write()

    if (channel == 'electron'): ch = '_el' 
    if (channel == 'muon'): ch = '_mu' 

    print channel,' data = ', VariablesSmooth['Data'+ch].Integral()
    print 'Total background = ', VariablesSmooth['TTbar_Madgraph'].Integral()+VariablesSmoothPre['WJets'].Integral()
    print 'WJets ',VariablesSmoothPre['WJets'].Integral()
    print 'TTbar ',VariablesSmooth['TTbar_Madgraph'].Integral()
    print 'Background / Data = ', (VariablesSmooth['TTbar_Madgraph'].Integral()+VariablesSmoothPre['WJets'].Integral())/VariablesSmooth['Data'+ch].Integral()
                
    ##### JES UP ##### 
    VariablesSmooth['WW_JESUP'].Add(VariablesSmooth['ZJets_M50_JESUP'])
    VariablesSmooth['WJets_JESUP'].Add(VariablesSmooth['WW_JESUP'])
    VariablesSmoothPre['WW_JESUP'].Add(VariablesSmoothPre['ZJets_M50_JESUP'])
    VariablesSmoothPre['WJets_JESUP'].Add(VariablesSmoothPre['WW_JESUP'])
    VariablesSmoothPre['WJets_JESUP'].SetName(chan+"invmass_"+btags+"__wjets__jes__plus")
    VariablesSmoothPre['WJets_JESUP'].SetTitle(chan+"invmass_"+btags+"__wjets__jes__plus")
    ######################################### 
    ### RESCALE W+Jets JES UP   
    #########################################
    print 'WJetsJESup ge1/pre: ',VariablesSmooth['WJets_JESUP'].Integral()/VariablesSmoothPre['WJets_JESUP'].Integral()
    VariablesSmoothPre['WJets_JESUP'].Scale( VariablesSmooth['WJets_JESUP'].Integral()/VariablesSmoothPre['WJets_JESUP'].Integral())
    VariablesSmoothPre['WJets_JESUP'].Write()

    VariablesSmooth['TTbar_Madgraph_JESUP'].SetName(chan+"invmass_"+btags+"__ttbar__jes__plus")
    VariablesSmooth['TTbar_Madgraph_JESUP'].SetTitle(chan+"invmass_"+btags+"__ttbar__jes__plus")
    VariablesSmooth['T_t_JESUP'].Add(VariablesSmooth['Tbar_t_JESUP'])
    VariablesSmooth['T_t_JESUP'].Add(VariablesSmooth['T_tW_JESUP'])
    VariablesSmooth['T_t_JESUP'].Add(VariablesSmooth['Tbar_tW_JESUP'])
    VariablesSmooth['TTbar_Madgraph_JESUP'].Add( VariablesSmooth['T_t_JESUP'])

    VariablesSmooth['T_s_JESUP'].Add(VariablesSmooth['Tbar_s_JESUP']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_JESUP'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_JESUP'].Write()
            VariablesSmooth['T_s_JESUP'].SetName(chan+"invmass_"+btags+"__topstb__jes__plus")
            VariablesSmooth['T_s_JESUP'].SetTitle(chan+"invmass_"+btags+"__topstb__jes__plus")
            VariablesSmooth['T_s_JESUP'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_JESUP'].Add(VariablesSmooth['T_s_JESUP'])
            VariablesSmooth['TTbar_Madgraph_JESUP'].Write()

    ##### JES DOWN #####
    VariablesSmooth['WW_JESDOWN'].Add(VariablesSmooth['ZJets_M50_JESDOWN'])
    VariablesSmooth['WJets_JESDOWN'].Add(VariablesSmooth['WW_JESDOWN'])
    VariablesSmoothPre['WW_JESDOWN'].Add(VariablesSmoothPre['ZJets_M50_JESDOWN'])
    VariablesSmoothPre['WJets_JESDOWN'].Add(VariablesSmoothPre['WW_JESDOWN'])
    VariablesSmoothPre['WJets_JESDOWN'].SetName(chan+"invmass_"+btags+"__wjets__jes__minus")
    VariablesSmoothPre['WJets_JESDOWN'].SetTitle(chan+"invmass_"+btags+"__wjets__jes__minus")
    ######################################### 
    ### REBIN W+Jets JES UP   
    #########################################
    print 'WJetsJESdown ge1/pre: ',VariablesSmooth['WJets_JESDOWN'].Integral()/VariablesSmoothPre['WJets_JESDOWN'].Integral()
    VariablesSmoothPre['WJets_JESDOWN'].Scale(VariablesSmooth['WJets_JESDOWN'].Integral()/VariablesSmoothPre['WJets_JESDOWN'].Integral())
    VariablesSmoothPre['WJets_JESDOWN'].Write()

    VariablesSmooth['TTbar_Madgraph_JESDOWN'].SetName(chan+"invmass_"+btags+"__ttbar__jes__minus")
    VariablesSmooth['TTbar_Madgraph_JESDOWN'].SetTitle(chan+"invmass_"+btags+"__ttbar__jes__minus")
    VariablesSmooth['T_t_JESDOWN'].Add(VariablesSmooth['Tbar_t_JESDOWN'])
    VariablesSmooth['T_t_JESDOWN'].Add(VariablesSmooth['T_tW_JESDOWN'])
    VariablesSmooth['T_t_JESDOWN'].Add(VariablesSmooth['Tbar_tW_JESDOWN'])
    VariablesSmooth['TTbar_Madgraph_JESDOWN'].Add( VariablesSmooth['T_t_JESDOWN'])

    VariablesSmooth['T_s_JESDOWN'].Add(VariablesSmooth['Tbar_s_JESDOWN']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_JESDOWN'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_JESDOWN'].Write()
            VariablesSmooth['T_s_JESDOWN'].SetName(chan+"invmass_"+btags+"__topstb__jes__minus")
            VariablesSmooth['T_s_JESDOWN'].SetTitle(chan+"invmass_"+btags+"__topstb__jes__minus")
            VariablesSmooth['T_s_JESDOWN'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_JESDOWN'].Add(VariablesSmooth['T_s_JESDOWN'])
            VariablesSmooth['TTbar_Madgraph_JESDOWN'].Write()

    ##### BTAG UP #####
    VariablesSmooth['WW_BTAGUP'].Add(VariablesSmooth['ZJets_M50_BTAGUP'])
    VariablesSmooth['WJets_BTAGUP'].Add(VariablesSmooth['WW_BTAGUP'])
    VariablesSmoothPre['WJets_BTAGUP'] = VariablesSmoothPre['WJets'].Clone()
    VariablesSmoothPre['WJets_BTAGUP'].SetName(chan+"invmass_"+btags+"__wjets__btag__plus")
    VariablesSmoothPre['WJets_BTAGUP'].SetTitle(chan+"invmass_"+btags+"__wjets__btag__plus")
    ######################################### 
    ### RESCALE W+Jets BTAG BASED ON RATIO
    #########################################
    VariablesSmoothPre['WJets_BTAGUP'].Scale(VariablesSmooth['WJets_BTAGUP'].Integral()/VariablesSmoothPre['WJets'].Integral() )
    VariablesSmoothPre['WJets_BTAGUP'].Write()

    VariablesSmooth['TTbar_Madgraph_BTAGUP'].SetName(chan+"invmass_"+btags+"__ttbar__btag__plus")
    VariablesSmooth['TTbar_Madgraph_BTAGUP'].SetTitle(chan+"invmass_"+btags+"__ttbar__btag__plus")
    VariablesSmooth['T_t_BTAGUP'].Add(VariablesSmooth['Tbar_t_BTAGUP'])
    VariablesSmooth['T_t_BTAGUP'].Add(VariablesSmooth['T_tW_BTAGUP'])
    VariablesSmooth['T_t_BTAGUP'].Add(VariablesSmooth['Tbar_tW_BTAGUP'])
    VariablesSmooth['TTbar_Madgraph_BTAGUP'].Add( VariablesSmooth['T_t_BTAGUP'])

    VariablesSmooth['T_s_BTAGUP'].Add(VariablesSmooth['Tbar_s_BTAGUP']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_BTAGUP'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_BTAGUP'].Write()
            VariablesSmooth['T_s_BTAGUP'].SetName(chan+"invmass_"+btags+"__topstb__btag__plus")
            VariablesSmooth['T_s_BTAGUP'].SetTitle(chan+"invmass_"+btags+"__topstb__btag__plus")
            VariablesSmooth['T_s_BTAGUP'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_BTAGUP'].Add(VariablesSmooth['T_s_BTAGUP'])
            VariablesSmooth['TTbar_Madgraph_BTAGUP'].Write()

    ##### BTAG DOWN #####
    VariablesSmooth['WW_BTAGDOWN'].Add(VariablesSmooth['ZJets_M50_BTAGDOWN'])
    VariablesSmooth['WJets_BTAGDOWN'].Add(VariablesSmooth['WW_BTAGDOWN'])
    VariablesSmoothPre['WJets_BTAGDOWN'] = VariablesSmoothPre['WJets'].Clone()
    VariablesSmoothPre['WJets_BTAGDOWN'].SetName(chan+"invmass_"+btags+"__wjets__btag__minus")
    VariablesSmoothPre['WJets_BTAGDOWN'].SetTitle(chan+"invmass_"+btags+"__wjets__btag__minus")
    ######################################### 
    ### RESCALE W+Jets BTAG BASED ON RATIO
    #########################################
    VariablesSmoothPre['WJets_BTAGDOWN'].Scale(VariablesSmooth['WJets_BTAGDOWN'].Integral()/VariablesSmoothPre['WJets'].Integral())
    VariablesSmoothPre['WJets_BTAGDOWN'].Write()

    VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].SetName(chan+"invmass_"+btags+"__ttbar__btag__minus")
    VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].SetTitle(chan+"invmass_"+btags+"__ttbar__btag__minus")
    VariablesSmooth['T_t_BTAGDOWN'].Add(VariablesSmooth['Tbar_t_BTAGDOWN'])
    VariablesSmooth['T_t_BTAGDOWN'].Add(VariablesSmooth['T_tW_BTAGDOWN'])
    VariablesSmooth['T_t_BTAGDOWN'].Add(VariablesSmooth['Tbar_tW_BTAGDOWN'])
    VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].Add( VariablesSmooth['T_t_BTAGDOWN'])

    VariablesSmooth['T_s_BTAGDOWN'].Add(VariablesSmooth['Tbar_s_BTAGDOWN']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        dummy = 1.0
        VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].Write()
            VariablesSmooth['T_s_BTAGDOWN'].SetName(chan+"invmass_"+btags+"__topstb__btag__minus")
            VariablesSmooth['T_s_BTAGDOWN'].SetTitle(chan+"invmass_"+btags+"__topstb__btag__minus")
            VariablesSmooth['T_s_BTAGDOWN'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].Add(VariablesSmooth['T_s_BTAGDOWN'])
            VariablesSmooth['TTbar_Madgraph_BTAGDOWN'].Write()

    ##### JER UP #####
    VariablesSmooth['WW_JERUP'].Add(VariablesSmooth['ZJets_M50_JERUP'])
    VariablesSmooth['WJets_JERUP'].Add(VariablesSmooth['WW_JERUP'])
    VariablesSmoothPre['WW_JERUP'].Add(VariablesSmoothPre['ZJets_M50_JERUP'])
    VariablesSmoothPre['WJets_JERUP'].Add(VariablesSmoothPre['WW_JERUP'])
    VariablesSmoothPre['WJets_JERUP'].SetName(chan+"invmass_"+btags+"__wjets__jer__plus")
    VariablesSmoothPre['WJets_JERUP'].SetTitle(chan+"invmass_"+btags+"__wjets__jer__plus")
    #########################################
    ### RESCALE W+Jets JER UP
    #########################################
    VariablesSmoothPre['WJets_JERUP'].Scale(VariablesSmooth['WJets_JERUP'].Integral()/VariablesSmoothPre['WJets_JERUP'].Integral())
    VariablesSmoothPre['WJets_JERUP'].Write()

    VariablesSmooth['TTbar_Madgraph_JERUP'].SetName(chan+"invmass_"+btags+"__ttbar__jer__plus")
    VariablesSmooth['TTbar_Madgraph_JERUP'].SetTitle(chan+"invmass_"+btags+"__ttbar__jer__plus")
    VariablesSmooth['T_t_JERUP'].Add(VariablesSmooth['Tbar_t_JERUP'])
    VariablesSmooth['T_t_JERUP'].Add(VariablesSmooth['T_tW_JERUP'])
    VariablesSmooth['T_t_JERUP'].Add(VariablesSmooth['Tbar_tW_JERUP'])
    VariablesSmooth['TTbar_Madgraph_JERUP'].Add( VariablesSmooth['T_t_JERUP'])

    VariablesSmooth['T_s_JERUP'].Add(VariablesSmooth['Tbar_s_JERUP']) 
   
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_JERUP'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_JERUP'].Write()
            VariablesSmooth['T_s_JERUP'].SetName(chan+"invmass_"+btags+"__topstb__jer__plus")
            VariablesSmooth['T_s_JERUP'].SetTitle(chan+"invmass_"+btags+"__topstb__jer__plus")
            VariablesSmooth['T_s_JERUP'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_JERUP'].Add(VariablesSmooth['T_s_JERUP'])
            VariablesSmooth['TTbar_Madgraph_JERUP'].Write()
 
    ##### JER DOWN #####
    VariablesSmooth['WW_JERDOWN'].Add(VariablesSmooth['ZJets_M50_JERDOWN'])
    VariablesSmooth['WJets_JERDOWN'].Add(VariablesSmooth['WW_JERDOWN'])
    VariablesSmoothPre['WW_JERDOWN'].Add(VariablesSmoothPre['ZJets_M50_JERDOWN'])
    VariablesSmoothPre['WJets_JERDOWN'].Add(VariablesSmoothPre['WW_JERDOWN'])
    VariablesSmoothPre['WJets_JERDOWN'].SetName(chan+"invmass_"+btags+"__wjets__jer__minus")
    VariablesSmoothPre['WJets_JERDOWN'].SetTitle(chan+"invmass_"+btags+"__wjets__jer__minus")
    #########################################
    ### RESCALE W+Jets JER DOWN   
    #########################################
    VariablesSmoothPre['WJets_JERDOWN'].Scale( VariablesSmooth['WJets_JERDOWN'].Integral()/VariablesSmoothPre['WJets_JERDOWN'].Integral())
    VariablesSmoothPre['WJets_JERDOWN'].Write()

    VariablesSmooth['TTbar_Madgraph_JERDOWN'].SetName(chan+"invmass_"+btags+"__ttbar__jer__minus")
    VariablesSmooth['TTbar_Madgraph_JERDOWN'].SetTitle(chan+"invmass_"+btags+"__ttbar__jer__minus")
    VariablesSmooth['T_t_JERDOWN'].Add(VariablesSmooth['Tbar_t_JERDOWN'])
    VariablesSmooth['T_t_JERDOWN'].Add(VariablesSmooth['T_tW_JERDOWN'])
    VariablesSmooth['T_t_JERDOWN'].Add(VariablesSmooth['Tbar_tW_JERDOWN'])
    VariablesSmooth['TTbar_Madgraph_JERDOWN'].Add( VariablesSmooth['T_t_JERDOWN'])

    VariablesSmooth['T_s_JERDOWN'].Add(VariablesSmooth['Tbar_s_JERDOWN']) 
   
    if (wprime == 'Left' or wprime == 'Mix'):
        VariablesSmooth['TTbar_Madgraph_JERDOWN'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_JERDOWN'].Write()
            VariablesSmooth['T_s_JERDOWN'].SetName(chan+"invmass_"+btags+"__topstb__jer__minus")
            VariablesSmooth['T_s_JERDOWN'].SetTitle(chan+"invmass_"+btags+"__topstb__jer__minus")
            VariablesSmooth['T_s_JERDOWN'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_JERDOWN'].Add(VariablesSmooth['T_s_JERDOWN'])
            VariablesSmooth['TTbar_Madgraph_JERDOWN'].Write()


    VariablesSmoothPUup['WW'].Add(VariablesSmoothPUup['ZJets_M50'])
    VariablesSmoothPUup['WJets'].Add(VariablesSmoothPUup['WW'])
    VariablesSmoothPUup['WJets'].SetName(chan+"invmass_"+btags+"__wjets__PU__plus")
    VariablesSmoothPUup['WJets'].SetTitle(chan+"invmass_"+btags+"__wjets__PU__plus")
    VariablesSmoothPUup['WJets'].Write()

    VariablesSmoothPUup['TTbar_Madgraph'].SetName(chan+"invmass_"+btags+"__ttbar__PU__plus")
    VariablesSmoothPUup['TTbar_Madgraph'].SetTitle(chan+"invmass_"+btags+"__ttbar__PU__plus")
    VariablesSmoothPUup['T_t'].Add(VariablesSmoothPUup['Tbar_t'])
    VariablesSmoothPUup['T_t'].Add(VariablesSmoothPUup['T_tW'])
    VariablesSmoothPUup['T_t'].Add(VariablesSmoothPUup['Tbar_tW'])
    VariablesSmoothPUup['TTbar_Madgraph'].Add( VariablesSmoothPUup['T_t'])

    VariablesSmoothPUup['T_s'].Add(VariablesSmoothPUup['Tbar_s']) 
   
    if (wprime == 'Left' or wprime == 'Mix'):
        VariablesSmoothPUup['TTbar_Madgraph'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmoothPUup['TTbar_Madgraph'].Write()
            VariablesSmoothPUup['T_s'].SetName(chan+"invmass_"+btags+"__topstb__PU__plus")
            VariablesSmoothPUup['T_s'].SetTitle(chan+"invmass_"+btags+"__topstb__PU__plus")
            VariablesSmoothPUup['T_s'].Write()
        elif (wprime == 'Right'):
            VariablesSmoothPUup['TTbar_Madgraph'].Add(VariablesSmoothPUup['T_s'])
            VariablesSmoothPUup['TTbar_Madgraph'].Write()

    VariablesSmoothPUup['WJetsDown'] = VariablesSmoothPre['WJets'].Clone()
    VariablesSmoothPUup['WJetsDown'].SetName(chan+"invmass_"+btags+"__wjets__PU__minus")
    VariablesSmoothPUup['WJetsDown'].SetTitle(chan+"invmass_"+btags+"__wjets__PU__minus")
    VariablesSmoothPUup['WJetsDown'].Write()
    VariablesSmoothPUup['TTbar_MadgraphDown'] = VariablesSmooth['TTbar_Madgraph'].Clone()
    VariablesSmoothPUup['TTbar_MadgraphDown'].SetName(chan+"invmass_"+btags+"__ttbar__PU__minus")
    VariablesSmoothPUup['TTbar_MadgraphDown'].SetTitle(chan+"invmass_"+btags+"__ttbar__PU__minus")
    VariablesSmoothPUup['TTbar_MadgraphDown'].Write()
    VariablesSmoothPUup['T_sDown'] = VariablesSmooth['T_s'].Clone()
    VariablesSmoothPUup['T_sDown'].SetName(chan+"invmass_"+btags+"__topstb__PU__minus")
    VariablesSmoothPUup['T_sDown'].SetTitle(chan+"invmass_"+btags+"__topstb__PU__minus")
    if (wprime == 'ModRight'): VariablesSmoothPUup['T_sDown'].Write()
            
 
    if (wprime == 'Right'):
        masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2600','2700','2800','2900','3000']
    else:
        masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2700','2800','2900']
    for mass in masses:
        if (wprime != 'ModRight'):
            VariablesSmooth['Wprime' + mass + wprime].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_JESUP'].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_JESDOWN'].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_JERUP'].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_JERDOWN'].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_BTAGUP'].Write()
            VariablesSmooth['Wprime' + mass + wprime + '_BTAGDOWN'].Write()
            VariablesSmoothPUup['Wprime' + mass + wprime ].Write()  
            VariablesSmoothPUup['Wprime' + mass + wprime+'Down' ] = VariablesSmooth['Wprime' + mass + wprime].Clone()
            VariablesSmoothPUup['Wprime' + mass + wprime+'Down' ].SetName( chan+'invmass_'+btags+'__'+'wp'+mass+'__PU__minus')
            VariablesSmoothPUup['Wprime' + mass + wprime+'Down' ].SetTitle( chan+'invmass_'+btags+'__'+'wp'+mass+'__PU__minus')
            VariablesSmoothPUup['Wprime' + mass + wprime+'Down' ].Write()
        else:
            VariablesSmooth['Wprime' + mass + 'Right'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_JESUP'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_JESDOWN'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_JERUP'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_JERDOWN'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_BTAGUP'].Write()
            VariablesSmooth['Wprime' + mass + 'Right' + '_BTAGDOWN'].Write()
            VariablesSmoothPUup['Wprime' + mass + 'Right' ].Write()  
            VariablesSmoothPUup['Wprime' + mass + 'Right'+'Down' ] = VariablesSmooth['Wprime' + mass + 'Right'].Clone()
            VariablesSmoothPUup['Wprime' + mass + 'Right'+'Down' ].SetName( chan+'invmass_'+btags+'__'+'wp'+mass+'__PU__minus')
            VariablesSmoothPUup['Wprime' + mass + 'Right'+'Down' ].SetTitle( chan+'invmass_'+btags+'__'+'wp'+mass+'__PU__minus')
            VariablesSmoothPUup['Wprime' + mass + 'Right'+'Down' ].Write() 

    
    VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].Add(VariablesSmooth['T_t'])
    VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].SetName(chan+"invmass_"+btags+"__ttbar__matching__minus")
    VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].SetTitle(chan+"invmass_"+btags+"__ttbar__matching__minus")
    if (wprime == 'Left' or wprime == 'Mix'):
        VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].Add(VariablesSmooth['T_s'])
            VariablesSmooth['TTbar_Madgraph_MATCHINGDOWN'].Write()

    VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].Add(VariablesSmooth['T_t'])
    VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].SetName(chan+"invmass_"+btags+"__ttbar__matching__plus")
    VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].SetTitle(chan+"invmass_"+btags+"__ttbar__matching__plus")
    if (wprime == 'Left' or wprime == 'Mix'):
        VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].Add(VariablesSmooth['T_s'])
            VariablesSmooth['TTbar_Madgraph_MATCHINGUP'].Write()

    VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].Add(VariablesSmooth['T_t'])
    VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].SetName(chan+"invmass_"+btags+"__ttbar__q2scale__minus")
    VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].SetTitle(chan+"invmass_"+btags+"__ttbar__q2scale__minus")
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].Add(VariablesSmooth['T_s'])
            VariablesSmooth['TTbar_Madgraph_SCALEDOWN'].Write()

    VariablesSmooth['TTbar_Madgraph_SCALEUP'].Add(VariablesSmooth['T_t'])
    VariablesSmooth['TTbar_Madgraph_SCALEUP'].SetName(chan+"invmass_"+btags+"__ttbar__q2scale__plus")
    VariablesSmooth['TTbar_Madgraph_SCALEUP'].SetTitle(chan+"invmass_"+btags+"__ttbar__q2scale__plus")
    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmooth['TTbar_Madgraph_SCALEUP'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmooth['TTbar_Madgraph_SCALEUP'].Write()
        elif (wprime == 'Right'):
            VariablesSmooth['TTbar_Madgraph_SCALEUP'].Add(VariablesSmooth['T_s'])
            VariablesSmooth['TTbar_Madgraph_SCALEUP'].Write()

    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['TTbar_Madgraph'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['T_t'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['Tbar_t'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['T_tW'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['Tbar_tW'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['T_s'],-1)
    VariablesSmooth0tag['Data'+ch].Add(VariablesSmooth0tag['Tbar_s'],-1)
    VariablesSmooth0tag['Data'+ch].Scale(VariablesSmooth['WJets'].Integral()/VariablesSmooth0tag['Data'+ch].Integral())
    VariablesSmooth0tag['Data'+ch].SetName(chan+"invmass_"+btags+"__wjets__zerotagshape__plus")
    VariablesSmooth0tag['Data'+ch].SetTitle(chan+"invmass_"+btags+"__wjets__zerotagshape__plus")
    for x in range(1,bin+1):  
        if (VariablesSmooth0tag['Data'+ch].GetBinContent(x) < 0): 
            VariablesSmooth0tag['Data'+ch].SetBinContent(x,0.000001)
        VariablesSmooth0tag['Data'+ch].SetBinError(x, math.sqrt(VariablesSmooth0tag['Data'+ch].GetBinContent(x)) )
    VariablesSmooth0tag['Data'+ch].Write()
    VariablesSmooth0tag['WJetsDown'] = VariablesSmoothPre['WJets'].Clone()
    VariablesSmooth0tag['WJetsDown'].SetName(chan+"invmass_"+btags+"__wjets__zerotagshape__minus")
    VariablesSmooth0tag['WJetsDown'].SetTitle(chan+"invmass_"+btags+"__wjets__zerotagshape__minus")
    VariablesSmooth0tag['WJetsDown'].Write()

    VariablesSmoothHFup['WW'].Add(VariablesSmoothHFup['ZJets_M50'])
    VariablesSmoothHFup['WJets'].Add(VariablesSmoothHFup['WW'])
    VariablesSmoothHFup['WJets'].SetName(chan+"invmass_"+btags+"__wjets__hf__plus")
    VariablesSmoothHFup['WJets'].SetTitle(chan+"invmass_"+btags+"__wjets__hf__plus")
    VariablesSmoothHFup['WJets'].Write()
    VariablesSmoothHFdown['WW'].Add(VariablesSmoothHFdown['ZJets_M50'])
    VariablesSmoothHFdown['WJets'].Add(VariablesSmoothHFdown['WW'])
    VariablesSmoothHFdown['WJets'].SetName(chan+"invmass_"+btags+"__wjets__hf__minus")
    VariablesSmoothHFdown['WJets'].SetTitle(chan+"invmass_"+btags+"__wjets__hf__minus")
    VariablesSmoothHFdown['WJets'].Write()

    VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].Add(VariablesSmooth['T_t'])
    VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].SetName(chan+"invmass_"+btags+"__ttbar__shape__plus")
    VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].SetTitle(chan+"invmass_"+btags+"__ttbar__shape__plus")
    VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].Add(VariablesSmooth['T_t'])
    VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].SetName(chan+"invmass_"+btags+"__ttbar__shape__minus")
    VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].SetTitle(chan+"invmass_"+btags+"__ttbar__shape__minus")

    if (wprime == 'Left' or wprime == 'Mix'): 
        VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].Write()
        VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].Write()
    else:
        if ( wprime == 'ModRight' ):
            VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].Write()
            VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].Write()
        elif (wprime == 'Right'):
            VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].Add(VariablesSmooth['T_s'])
            VariablesSmoothTTbarShapeUp['TTbar_Madgraph'].Write()
            VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].Add(VariablesSmooth['T_s'])
            VariablesSmoothTTbarShapeDown['TTbar_Madgraph'].Write()



  
#wprime = 'Right'
var = 'BestJetJet2W_M_LjetsTopoCalcNew'; high = 4000; xaxis = "W' invariant mass [GeV/c^{2}]"; yaxis = 'Events / 10 GeV'; save = 'BestJetJet2W_M'

#btags = 0
#btagstr = '0btags'
btags = 'final'
xlow  = [ 100, 300,  400, 500, 600, 700, 800, 900,  1000, 1100, 1200,  1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 4000]
#xlow  = [ 100, 300,  400, 500, 600, 700, 800, 900,  1000, 1100, 1200,  1300, 1400, 1500, 1600, 1700, 1800, 4000]
#xlow  = [ 100, 300,  400, 500, 600, 700, 800, 900,  1000, 1100, 1200,  1300, 1400, 1500, 1600, 1700, 1800, 2000, 4000]


bins = len(xlow)-1 

List_DataBgEl_RightWprime = copy.copy(List_DataEl) 
List_DataBgEl_RightWprime.extend(List_Bg) 
List_DataBgEl_RightWprime.extend(List_Right) 

List_DataBgEl_ModRightWprime = copy.copy(List_DataEl) 
List_DataBgEl_ModRightWprime.extend(List_Bg) 
List_DataBgEl_ModRightWprime.extend(List_ModRight) 

List_DataBgEl_LeftWprime = copy.copy(List_DataEl) 
List_DataBgEl_LeftWprime.extend(List_Bg) 
List_DataBgEl_LeftWprime.extend(List_Left) 

List_DataBgEl_MixWprime = copy.copy(List_DataEl) 
List_DataBgEl_MixWprime.extend(List_Bg) 
List_DataBgEl_MixWprime.extend(List_Mix) 

List_DataBgMu_RightWprime = copy.copy(List_DataMu) 
List_DataBgMu_RightWprime.extend(List_Bg) 
List_DataBgMu_RightWprime.extend(List_Right) 

List_DataBgMu_ModRightWprime = copy.copy(List_DataMu) 
List_DataBgMu_ModRightWprime.extend(List_Bg) 
List_DataBgMu_ModRightWprime.extend(List_ModRight) 

List_DataBgMu_LeftWprime = copy.copy(List_DataMu) 
List_DataBgMu_LeftWprime.extend(List_Bg) 
List_DataBgMu_LeftWprime.extend(List_Left) 

List_DataBgMu_MixWprime = copy.copy(List_DataMu) 
List_DataBgMu_MixWprime.extend(List_Bg) 
List_DataBgMu_MixWprime.extend(List_Mix) 


channel = 'electron'

wprime = 'ModRight'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_RightWprime) 

wprime = 'Right'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_RightWprime) 

wprime = 'Left'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_LeftWprime) 

wprime = 'Mix'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgEl_MixWprime) 


channel = 'muon'

wprime = 'ModRight'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_RightWprime) 

wprime = 'Right'
f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
f.cd() 
plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_RightWprime) 

wprime = 'Left'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_LeftWprime) 

wprime = 'Mix'
#f = TFile("RootFiles_For2DLimits_06Mar_finalbins_scaleTTbarTopPtEMU750_landau/"+channel+"_"+save+"_Wprime"+wprime+"_Histos-"+btags+"_85_140_dr03_lep50.root","RECREATE")
#f.cd() 
#plot_DataVsMc(channel,var, '85', '140', '0.3', bins, xlow, high, yaxis, xaxis , save, wprime, btags, List_DataBgMu_MixWprime) 

