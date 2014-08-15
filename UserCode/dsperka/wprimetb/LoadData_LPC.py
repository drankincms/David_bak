from ROOT import TFile

#lumi_el = 19624.0 
#lumi_mu = 19624.0 
lumi_el = 19521.0
lumi_mu = 19521.0

lumiPlot_el= '19.5'
lumiPlot_mu= '19.5'

#outdir = "/uscms_data/d2/dsperka/8TeV/Samples/07Feb_All/"
outdir = "/uscms_data/d2/dsperka/8TeV/Samples/29May_All/"
outdir2 = "/uscms_data/d2/dsperka/8TeV/Samples/29May_All/"

######################################################
RootFiles = {}
Trees = {}

RootFiles['Data_el'] = TFile(outdir+"Data_el_19pt6fb.root") # called 19pt6, but its 19pt5 with pixelLumi 
RootFiles['Data_mu'] = TFile(outdir+"Data_mu_19pt6fb.root") # called 19pt6, but its 19pt5 with pixelLumi
Trees['Data_el']  = RootFiles['Data_el'].Get("ljmet")
Trees['Data_mu']  = RootFiles['Data_mu'].Get("ljmet")

Backgrounds = ['WJets','W1Jets','W2Jets','W3Jets','W4Jets','WW','TTbar_Madgraph','TTbar_Powheg','ZJets_M50','T_t','Tbar_t','T_tW','Tbar_tW','T_s','Tbar_s','QCD_Pt_80_170_EM','QCD_Pt_170_250_EM','QCD_Pt_250_350_EM','QCD_Pt_350_EM']
for bg in Backgrounds:
    RootFiles[bg] = TFile(outdir2+bg+".root")
    RootFiles[bg+'_JESUP'] = TFile(outdir2+"JESUP/"+bg+"_JESUP.root")
    RootFiles[bg+'_JESDOWN'] = TFile(outdir2+"JESDOWN/"+bg+"_JESDOWN.root")
    RootFiles[bg+'_JERUP'] = TFile(outdir2+"JERUP/"+bg+"_JERUP.root")
    RootFiles[bg+'_JERDOWN'] = TFile(outdir2+"JERDOWN/"+bg+"_JERDOWN.root")
    RootFiles[bg+'_BTAGUP'] = TFile(outdir2+"BTAGUP/"+bg+"_BTAGUP.root")
    RootFiles[bg+'_BTAGDOWN'] = TFile(outdir2+"BTAGDOWN/"+bg+"_BTAGDOWN.root")

    Trees[bg] = RootFiles[bg].Get("ljmet")
    Trees[bg+'_JESUP'] = RootFiles[bg+'_JESUP'].Get("ljmet")
    Trees[bg+'_JESDOWN'] = RootFiles[bg+'_JESDOWN'].Get("ljmet")
    Trees[bg+'_JERUP'] = RootFiles[bg+'_JERUP'].Get("ljmet")
    Trees[bg+'_JERDOWN'] = RootFiles[bg+'_JERDOWN'].Get("ljmet")
    Trees[bg+'_BTAGUP'] = RootFiles[bg+'_BTAGUP'].Get("ljmet")
    Trees[bg+'_BTAGDOWN'] = RootFiles[bg+'_BTAGDOWN'].Get("ljmet")

RootFiles['TTbar_Madgraph_MATCHINGUP'] = TFile(outdir2+"TTbar_matchingup.root")
RootFiles['TTbar_Madgraph_MATCHINGDOWN'] = TFile(outdir2+"TTbar_matchingdown.root")
RootFiles['TTbar_Madgraph_SCALEUP'] = TFile(outdir2+"TTbar_scaleup.root")
RootFiles['TTbar_Madgraph_SCALEDOWN'] = TFile(outdir2+"TTbar_scaledown.root")

Trees['TTbar_Madgraph_MATCHINGUP'] = RootFiles['TTbar_Madgraph_MATCHINGUP'].Get("ljmet")
Trees['TTbar_Madgraph_MATCHINGDOWN'] = RootFiles['TTbar_Madgraph_MATCHINGDOWN'].Get("ljmet")
Trees['TTbar_Madgraph_SCALEUP'] = RootFiles['TTbar_Madgraph_SCALEUP'].Get("ljmet")
Trees['TTbar_Madgraph_SCALEDOWN'] = RootFiles['TTbar_Madgraph_SCALEDOWN'].Get("ljmet")

masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2600','2700','2800','2900','3000']

for mass in masses:
    for coup in ['Right','Left','Mix']:
        RootFiles['Wprime'+mass+coup] = TFile(outdir+"Wprime"+mass+coup+".root")
        RootFiles['Wprime'+mass+coup+'_JESUP'] = TFile(outdir+"JESUP/Wprime"+mass+coup+"_JESUP.root")
        RootFiles['Wprime'+mass+coup+'_JESDOWN'] = TFile(outdir+"JESDOWN/Wprime"+mass+coup+"_JESDOWN.root")
        RootFiles['Wprime'+mass+coup+'_JERUP'] = TFile(outdir+"JERUP/Wprime"+mass+coup+"_JERUP.root")
        RootFiles['Wprime'+mass+coup+'_JERDOWN'] = TFile(outdir+"JERDOWN/Wprime"+mass+coup+"_JERDOWN.root")
        RootFiles['Wprime'+mass+coup+'_BTAGUP'] = TFile(outdir+"BTAGUP/Wprime"+mass+coup+"_BTAGUP.root")
        RootFiles['Wprime'+mass+coup+'_BTAGDOWN'] = TFile(outdir+"BTAGDOWN/Wprime"+mass+coup+"_BTAGDOWN.root")

        Trees['Wprime'+mass+coup] = RootFiles['Wprime'+mass+coup].Get("ljmet")
        Trees['Wprime'+mass+coup+'_JESUP'] = RootFiles['Wprime'+mass+coup+'_JESUP'].Get("ljmet")
        Trees['Wprime'+mass+coup+'_JESDOWN'] = RootFiles['Wprime'+mass+coup+'_JESDOWN'].Get("ljmet")
        Trees['Wprime'+mass+coup+'_JERUP'] = RootFiles['Wprime'+mass+coup+'_JERUP'].Get("ljmet")
        Trees['Wprime'+mass+coup+'_JERDOWN'] = RootFiles['Wprime'+mass+coup+'_JERDOWN'].Get("ljmet")
        Trees['Wprime'+mass+coup+'_BTAGUP'] = RootFiles['Wprime'+mass+coup+'_BTAGUP'].Get("ljmet")
        Trees['Wprime'+mass+coup+'_BTAGDOWN'] = RootFiles['Wprime'+mass+coup+'_BTAGDOWN'].Get("ljmet")


#############################################
Nevents = {}
Nevents['WJets'] = 76041475.0 #v1+v2
Nevents['W1Jets'] = 23141598.0
Nevents['W2Jets'] = 34044921.0
Nevents['W3Jets'] = 15539503.0
Nevents['W4Jets'] = 13382803.0
Nevents['WW'] = 10000431.0
Nevents['TTbar_Madgraph'] = 6923750.0 
Nevents['TTbar_Powheg'] = 21591169.0 
Nevents['TTbar_Total'] = (21591169.0+6923750.0) 
Nevents['ZJets_M50'] = 30459503.0
Nevents['T_tW'] = 497658.0
Nevents['Tbar_tW'] = 493460.0
Nevents['T_t'] = 3758227.0
Nevents['Tbar_t'] = 1935072.0
Nevents['T_s'] = 259961.0
Nevents['Tbar_s'] = 139974.0
Nevents['QCD_Pt_80_170_EM'] =  34542763.0 
Nevents['QCD_Pt_170_250_EM'] = 31697066.0
Nevents['QCD_Pt_250_350_EM'] = 34611322.0
Nevents['QCD_Pt_350_EM'] = 34080562.0


Nevents['TTbar_Madgraph_MATCHINGUP'] = 5415010.0
Nevents['TTbar_Madgraph_MATCHINGDOWN'] = 5476728.0
Nevents['TTbar_Madgraph_SCALEUP'] = 5009488.0
Nevents['TTbar_Madgraph_SCALEDOWN'] = 5387181.0

Nevents['Wprime800Right'] = 920654.0
Nevents['Wprime900Right'] = 942816.0
Nevents['Wprime1000Right'] = 907958.0
Nevents['Wprime1100Right'] = 831508.0
Nevents['Wprime1200Right'] = 965528.0
Nevents['Wprime1300Right'] = 881046.0
Nevents['Wprime1400Right'] = 920262.0
Nevents['Wprime1500Right'] = 907297.0
Nevents['Wprime1600Right'] = 892146.0
Nevents['Wprime1700Right'] = 924438.0
Nevents['Wprime1800Right'] = 841448.0
Nevents['Wprime1900Right'] = 835381.0
Nevents['Wprime2000Right'] = 841836.0
Nevents['Wprime2100Right'] = 926108.0
Nevents['Wprime2200Right'] = 932785.0
Nevents['Wprime2300Right'] = 784768.0
Nevents['Wprime2400Right'] = 894786.0
Nevents['Wprime2500Right'] = 878643.0
Nevents['Wprime2600Right'] = 944599.0
Nevents['Wprime2700Right'] = 915158.0
Nevents['Wprime2800Right'] = 835281.0
Nevents['Wprime2900Right'] = 910111.0
Nevents['Wprime3000Right'] = 932601.0

Nevents['Wprime800Left'] = 941306.0
Nevents['Wprime900Left'] = 906657.0
Nevents['Wprime1000Left'] = 908337.0
Nevents['Wprime1100Left'] = 798919.0
Nevents['Wprime1200Left'] = 959534.0
Nevents['Wprime1300Left'] = 963820.0
Nevents['Wprime1400Left'] = 942066.0
Nevents['Wprime1500Left'] = 952749.0
Nevents['Wprime1600Left'] = 954829.0
Nevents['Wprime1700Left'] = 948063.0
Nevents['Wprime1800Left'] = 936673.0
Nevents['Wprime1900Left'] = 911699.0
#Nevents['Wprime2000Left'] = 926534.0
Nevents['Wprime2000Left'] = 903371.0
#Nevents['Wprime2100Left'] = 931047.0
Nevents['Wprime2100Left'] = 861801.0
#Nevents['Wprime2200Left'] = 945666.0
Nevents['Wprime2200Left'] = 922413.0
#Nevents['Wprime2300Left'] = 971884.0
Nevents['Wprime2300Left'] = 964171.0
Nevents['Wprime2400Left'] = 931031.0
Nevents['Wprime2500Left'] = 911826.0
Nevents['Wprime2600Left'] = 931038.0
Nevents['Wprime2700Left'] = 907930.0
Nevents['Wprime2800Left'] = 917514.0
Nevents['Wprime2900Left'] = 940379.0
Nevents['Wprime3000Left'] = 934903.0

Nevents['Wprime800Mix'] = 920851.0 
Nevents['Wprime900Mix'] = 962105.0 
Nevents['Wprime1000Mix'] = 952695.0 
Nevents['Wprime1100Mix'] = 499057.0
Nevents['Wprime1200Mix'] = 949408.0
Nevents['Wprime1300Mix'] = 957707.0
Nevents['Wprime1400Mix'] = 499049.0
Nevents['Wprime1500Mix'] = 972899.0
Nevents['Wprime1600Mix'] = 948242.0
Nevents['Wprime1700Mix'] = 951497.0
Nevents['Wprime1800Mix'] = 963803.0
Nevents['Wprime1900Mix'] = 978267.0
Nevents['Wprime2000Mix'] = 929173.0
Nevents['Wprime2100Mix'] = 913931.0
Nevents['Wprime2200Mix'] = 938946.0
Nevents['Wprime2300Mix'] = 903118.0
Nevents['Wprime2400Mix'] = 956188.0
Nevents['Wprime2500Mix'] = 962673.0
Nevents['Wprime2600Mix'] = 945159.0
Nevents['Wprime2700Mix'] = 919176.0
Nevents['Wprime2800Mix'] = 921391.0
Nevents['Wprime2900Mix'] = 927989.0
Nevents['Wprime3000Mix'] = 932353.0

for sys in ['JESUP','JESDOWN','JERUP','JERDOWN','BTAGUP','BTAGDOWN']:
    for bg in Backgrounds:
        Nevents[bg+'_'+sys] = Nevents[bg]
    for mass in masses:
        for coup in ['Right','Left','Mix']:
            print 'Wprime'+mass+coup+'_'+sys
            Nevents['Wprime'+mass+coup+'_'+sys] = Nevents['Wprime'+mass+coup]

#############################################
xsec = {}
xsec['WJets'] = 37503.0 # fixme, adding zjets number
xsec['W1Jets'] = 5400.0*1.23
xsec['W2Jets'] = 1750.0*1.23
xsec['W3Jets'] = 519.0*1.23
xsec['W4Jets'] = 214.0*1.23
xsec['WW'] = 54.838 # CTEQ
#xsec['TTbar_Madgraph'] = 225.197 # MCFM NLO
#xsec['TTbar_Powheg'] = 225.197 # MCFM NLO
#xsec['TTbar_Madgraph'] = 234.0 # approx NNLO
#xsec['TTbar_Powheg'] = 234.0 # approx NNLO
xsec['TTbar_Madgraph'] = 245.0 # NNLO
xsec['TTbar_Powheg'] = 245.0 # NNLO
xsec['TTbar_Total'] = 245.0 # NNLO
xsec['ZJets_M50'] = 3503.71
xsec['T_tW'] = 11.1
xsec['Tbar_tW'] = 11.1
xsec['T_t'] = 56.4
xsec['Tbar_t'] = 30.7
xsec['T_s'] = 3.79
xsec['Tbar_s'] = 1.76

xsec['TTbar_Madgraph_MATCHINGUP'] = 245.0
xsec['TTbar_Madgraph_MATCHINGDOWN'] = 245.0
xsec['TTbar_Madgraph_SCALEUP'] = 245.0
xsec['TTbar_Madgraph_SCALEDOWN'] = 245.0

xsec['QCD_Pt_80_170_EM'] = 1191000.0*0.1539
xsec['QCD_Pt_170_250_EM'] = 30990.0*0.148
xsec['QCD_Pt_250_350_EM'] = 4250.0*0.131
xsec['QCD_Pt_350_EM'] = 810.0*0.11

xsec['Wprime800Right'] = 2.3022
xsec['Wprime900Right'] = 1.3818
xsec['Wprime1000Right'] = 0.85538
xsec['Wprime1100Right'] = 0.54325
xsec['Wprime1200Right'] = 0.35203
xsec['Wprime1300Right'] = 0.23219
xsec['Wprime1400Right'] = 0.15547
xsec['Wprime1500Right'] = 0.10518
xsec['Wprime1600Right'] = 0.072012
xsec['Wprime1700Right'] = 0.049683
xsec['Wprime1800Right'] = 0.034576
xsec['Wprime1900Right'] = 0.024249
xsec['Wprime2000Right'] = 0.017124
xsec['Wprime2100Right'] = 0.012176
xsec['Wprime2200Right'] = 0.0087191
xsec['Wprime2300Right'] = 0.0062918
xsec['Wprime2400Right'] = 0.0045757
xsec['Wprime2500Right'] = 0.0033568
xsec['Wprime2600Right'] = 0.0024870
xsec['Wprime2700Right'] = 0.0018624
xsec['Wprime2800Right'] = 0.0014102
xsec['Wprime2900Right'] = 0.0010818
xsec['Wprime3000Right'] = 0.00084115

xsec['Wprime800Left'] = 3.1080
xsec['Wprime900Left'] = 2.2731
xsec['Wprime1000Left'] = 1.8087
xsec['Wprime1100Left'] = 1.547
xsec['Wprime1200Left'] = 1.3870
xsec['Wprime1300Left'] = 1.2945
xsec['Wprime1400Left'] = 1.2390
xsec['Wprime1500Left'] = 1.2061
xsec['Wprime1600Left'] = 1.1869
xsec['Wprime1700Left'] = 1.1761
xsec['Wprime1800Left'] = 1.1705
xsec['Wprime1900Left'] = 1.1678
xsec['Wprime2000Left'] = 1.1673
xsec['Wprime2100Left'] = 1.1680
xsec['Wprime2200Left'] = 1.1692
xsec['Wprime2300Left'] = 1.1711
xsec['Wprime2400Left'] = 1.1727
xsec['Wprime2500Left'] = 1.1746
xsec['Wprime2600Left'] = 1.1763
xsec['Wprime2700Left'] = 1.1780
xsec['Wprime2800Left'] = 1.1797
xsec['Wprime2900Left'] = 1.1810
xsec['Wprime3000Left'] = 1.1825

xsec['Wprime800Mix'] = 5.4166
xsec['Wprime900Mix'] = 3.6684
xsec['Wprime1000Mix'] = 2.6815
xsec['Wprime1100Mix'] = 2.1031
xsec['Wprime1200Mix'] = 1.7539
xsec['Wprime1300Mix'] = 1.5389
xsec['Wprime1400Mix'] = 1.4043
xsec['Wprime1500Mix'] = 1.3194
xsec['Wprime1600Mix'] = 1.2650
xsec['Wprime1700Mix'] = 1.2305
xsec['Wprime1800Mix'] = 1.2090
xsec['Wprime1900Mix'] = 1.1954
xsec['Wprime2000Mix'] = 1.1872
xsec['Wprime2100Mix'] = 1.1824
xsec['Wprime2200Mix'] = 1.1798
xsec['Wprime2300Mix'] = 1.1787
xsec['Wprime2400Mix'] = 1.1784
xsec['Wprime2500Mix'] = 1.1791
xsec['Wprime2600Mix'] = 1.1792
xsec['Wprime2700Mix'] = 1.1803
xsec['Wprime2800Mix'] = 1.1813
xsec['Wprime2900Mix'] = 1.1825
xsec['Wprime3000Mix'] = 1.1835

xsec_norm = {}

for bg in Backgrounds:
    for sys in ['JESUP','JESDOWN','JERUP','JERDOWN','BTAGUP','BTAGDOWN']:
        xsec[bg+'_'+sys] = xsec[bg]
        xsec_norm[bg] = xsec[bg]
        xsec_norm[bg+'_'+sys] = xsec[bg]

for mass in masses:
    for coup in ['Right','Left','Mix']:
        xsec_norm['Wprime'+mass+coup] = 1.0
        for sys in ['JESUP','JESDOWN','JERUP','JERDOWN','BTAGUP','BTAGDOWN']:
            xsec['Wprime'+mass+coup+'_'+sys] = xsec['Wprime'+mass+coup]
            xsec_norm['Wprime'+mass+coup+'_'+sys] = 1.0

xsec_norm['TTbar_Madgraph_MATCHINGUP'] = 245.0
xsec_norm['TTbar_Madgraph_MATCHINGDOWN'] = 245.0
xsec_norm['TTbar_Madgraph_SCALEUP'] = 245.0
xsec_norm['TTbar_Madgraph_SCALEDOWN'] = 245.0

#############################################
############### FOR BDT #####################
#############################################


BDTDir = 'BDTDone/'
chanlist = ['el','mu']
masslist = ['wp1700','wp1900','wp2100']
couplist = ['R']
Typelist = ['data','wjets','ttbar','t','bt','tw','btw','s','bs','zjets','ww']
taglist = ['GE1BTag']
Samples = []

############# Root Files ###################
RootFilesBDT = {}
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
            for Type in Typelist:
                #RootFilesBDT[mass+coup+'_'+Type+'_'+chan] = TFile(BDTDir+'BDT_GE1BTag_'+mass+coup+'_'+Type+'_'+chan+'.root')
                #print 'Loaded Root file for ',mass+coup+'_'+Type+'_'+chan,': ',BDTDir+'BDT_GE1BTag_'+mass+coup+'_'+Type+'_'+chan+'.root'
                Samples.extend([mass+coup+'_'+Type+'_'+chan])
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
            #RootFilesBDT[mass+coup+'_'+mass+coup+'_'+chan] = TFile(BDTDir+'BDT_GE1BTag_'+mass+coup+'_'+mass+coup+'_'+chan+'.root')
            #print 'Loaded Root file for ',mass+coup+'_'+mass+coup+'_'+chan,': ',BDTDir+'BDT_GE1BTag_'+mass+coup+'_'+mass+coup+'_'+chan+'.root'
            Samples.extend([mass+coup+'_'+mass+coup+'_'+chan])


############## Trees #######################
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
            for Type in Typelist:
                #Trees[mass+coup+'_'+Type+'_'+chan]  = RootFilesBDT[mass+coup+'_'+Type+'_'+chan].Get(Type)
                #print 'Loaded Tree for ',mass+coup+'_'+Type+'_'+chan,': ',mass+coup+'_'+Type+'_'+chan
                dummy=1.0
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
            #Trees[mass+coup+'_'+mass+coup+'_'+chan]  = RootFilesBDT[mass+coup+'_'+mass+coup+'_'+chan].Get(mass+coup)
            #print 'Loaded Tree for ',mass+coup+'_'+mass+coup+'_'+chan,': ',mass+coup+'_'+mass+coup+'_'+chan
            dummy=1.0


#### Generated Events (same as above) ######
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
                Nevents[mass+coup+'_wjets_'+chan] = Nevents['WJets']
                Nevents[mass+coup+'_ww_'+chan] = Nevents['WW']
                Nevents[mass+coup+'_ttbar_'+chan] = Nevents['TTbar_Madgraph']
                Nevents[mass+coup+'_zjets_'+chan] = Nevents['ZJets_M50']
                Nevents[mass+coup+'_tw_'+chan] = Nevents['T_tW']
                Nevents[mass+coup+'_btw_'+chan] = Nevents['Tbar_tW']
                Nevents[mass+coup+'_t_'+chan] = Nevents['T_t']
                Nevents[mass+coup+'_bt_'+chan] = Nevents['Tbar_t']
                Nevents[mass+coup+'_s_'+chan] = Nevents['T_s']
                Nevents[mass+coup+'_bs_'+chan] = Nevents['Tbar_s']

for chan in chanlist:
    Nevents['wp800R_wp800R_'+chan] = Nevents['Wprime800Right'] 
    Nevents['wp900R_wp900R_'+chan] = Nevents['Wprime900Right'] 
    Nevents['wp1000R_wp1000R_'+chan] = Nevents['Wprime1000Right'] 
    Nevents['wp1100R_wp1100R_'+chan] = Nevents['Wprime1100Right'] 
    Nevents['wp1200R_wp1200R_'+chan] = Nevents['Wprime1200Right'] 
    Nevents['wp1300R_wp1300R_'+chan] = Nevents['Wprime1300Right'] 
    Nevents['wp1400R_wp1400R_'+chan] = Nevents['Wprime1400Right'] 
    Nevents['wp1500R_wp1500R_'+chan] = Nevents['Wprime1500Right'] 
    Nevents['wp1600R_wp1600R_'+chan] = Nevents['Wprime1600Right'] 
    Nevents['wp1700R_wp1700R_'+chan] = Nevents['Wprime1700Right'] 
    Nevents['wp1800R_wp1800R_'+chan] = Nevents['Wprime1800Right'] 
    Nevents['wp1900R_wp1900R_'+chan] = Nevents['Wprime1900Right'] 
    Nevents['wp2000R_wp2000R_'+chan] = Nevents['Wprime2000Right'] 
    Nevents['wp2100R_wp2100R_'+chan] = Nevents['Wprime2100Right'] 
    Nevents['wp2200R_wp2200R_'+chan] = Nevents['Wprime2200Right'] 
    Nevents['wp2300R_wp2300R_'+chan] = Nevents['Wprime2300Right'] 
    Nevents['wp2400R_wp2400R_'+chan] = Nevents['Wprime2400Right'] 
    Nevents['wp2500R_wp2500R_'+chan] = Nevents['Wprime2500Right'] 
    Nevents['wp2600R_wp2600R_'+chan] = Nevents['Wprime2600Right'] 
    Nevents['wp2700R_wp2700R_'+chan] = Nevents['Wprime2700Right'] 
    Nevents['wp2800R_wp2800R_'+chan] = Nevents['Wprime2800Right'] 
    Nevents['wp2900R_wp2900R_'+chan] = Nevents['Wprime2900Right'] 
    Nevents['wp3000R_wp3000R_'+chan] = Nevents['Wprime3000Right'] 

#### Cross sections (same as above) ######
for chan in chanlist:
    for mass in masslist:
        for coup in couplist:
                xsec[mass+coup+'_wjets_'+chan] = xsec['WJets']
                xsec[mass+coup+'_ww_'+chan] = xsec['WW']
                xsec[mass+coup+'_ttbar_'+chan] = xsec['TTbar_Madgraph']
                xsec[mass+coup+'_zjets_'+chan] = xsec['ZJets_M50']
                xsec[mass+coup+'_tw_'+chan] = xsec['T_tW']
                xsec[mass+coup+'_btw_'+chan] = xsec['Tbar_tW']
                xsec[mass+coup+'_t_'+chan] = xsec['T_t']
                xsec[mass+coup+'_bt_'+chan] = xsec['Tbar_t']
                xsec[mass+coup+'_s_'+chan] = xsec['T_s']
                xsec[mass+coup+'_bs_'+chan] = xsec['Tbar_s']

for chan in chanlist:
    xsec['wp800R_wp800R_'+chan] = xsec['Wprime800Right'] 
    xsec['wp900R_wp900R_'+chan] = xsec['Wprime900Right'] 
    xsec['wp1000R_wp1000R_'+chan] = xsec['Wprime1000Right'] 
    xsec['wp1100R_wp1100R_'+chan] = xsec['Wprime1100Right'] 
    xsec['wp1200R_wp1200R_'+chan] = xsec['Wprime1200Right'] 
    xsec['wp1300R_wp1300R_'+chan] = xsec['Wprime1300Right'] 
    xsec['wp1400R_wp1400R_'+chan] = xsec['Wprime1400Right'] 
    xsec['wp1500R_wp1500R_'+chan] = xsec['Wprime1500Right'] 
    xsec['wp1600R_wp1600R_'+chan] = xsec['Wprime1600Right'] 
    xsec['wp1700R_wp1700R_'+chan] = xsec['Wprime1700Right'] 
    xsec['wp1800R_wp1800R_'+chan] = xsec['Wprime1800Right'] 
    xsec['wp1900R_wp1900R_'+chan] = xsec['Wprime1900Right'] 
    xsec['wp2000R_wp2000R_'+chan] = xsec['Wprime2000Right'] 
    xsec['wp2100R_wp2100R_'+chan] = xsec['Wprime2100Right'] 
    xsec['wp2200R_wp2200R_'+chan] = xsec['Wprime2200Right'] 
    xsec['wp2300R_wp2300R_'+chan] = xsec['Wprime2300Right'] 
    xsec['wp2400R_wp2400R_'+chan] = xsec['Wprime2400Right'] 
    xsec['wp2500R_wp2500R_'+chan] = xsec['Wprime2500Right'] 
    xsec['wp2600R_wp2600R_'+chan] = xsec['Wprime2600Right'] 
    xsec['wp2700R_wp2700R_'+chan] = xsec['Wprime2700Right'] 
    xsec['wp2800R_wp2800R_'+chan] = xsec['Wprime2800Right'] 
    xsec['wp2900R_wp2900R_'+chan] = xsec['Wprime2900Right'] 
    xsec['wp3000R_wp3000R_'+chan] = xsec['Wprime3000Right'] 


