#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 50

rel_base = os.environ['CMSSW_BASE']

outdir = '/uscms_data/d2/dsperka/8TeV/Samples/29May/'

### What is the name of your FWLite Analyzer
FWLiteAnalyzer = 'ljmet'

### Selection to run
SELECTOR = "'WprimeSelector'"

### Is this sample MC?
ISMCSAMPLE = 'True'

### Is this a WJets sample?
ISWJETS = 'False'

### Do 53X JEC?
DO53XJEC = 'True'

### Use best top for neutrino pz solution?
USEBESTTOP = 'True'
 
### Triggers
MCTRIGGEREL = "'HLT_Ele27_WP80_v10'"
MCTRIGGERMU = "'HLT_IsoMu24_eta2p1_v13'"
DATATRIGGEREL = "'HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'"
DATATRIGGERMU = "'HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'"

### Lepton Selection
MINTIGHTMUON = '0'
MINTIGHTELECTRON = '0'
MINTIGHTLEPTON = '1'
MAXTIGHTLEPTON = '1'

### Jet Selection
LEADINGJETPT = '80.0'

### Which Systematics to do
DONOMINAL = 'True'
DOJES = 'True'
DOJER = 'True'
DOBTAG = 'True'
DOQCDMC = 'True'
DOTTBARSYS = 'True'

### JSON file to use
MYJSON = "''"

### for laser calibration correction filter ###
DOLASERCALFILT = 'False'

### Systematics flags
BTAGUNCERTUP = 'False'
BTAGUNCERTDOWN = 'False'
JECUNCERTUP = 'False'
JECUNCERTDOWN = 'False'
JERUNCERTUP = 'False'
JERUNCERTDOWN = 'False'


#################################################
### Names to give to your output root files
#################################################

prefix = []

if DONOMINAL=='True':
    prefix.extend([
    'WJets',
    'W1Jets',
    'W2Jets',
    'W3Jets',
    'W4Jets',
    'TTbar_Madgraph',
    'TTbar_Powheg',
    'T_t',
    'Tbar_t',
    'T_tW',
    'Tbar_tW',
    'T_s',
    'Tbar_s',
    'WW',
    'ZJets_M50',
    'ZZ',
    #'TT_Mtt_700to1000',
    #'TT_Mtt_1000toInf',
    ])

    if DOJES=='True':
        prefix.extend([
    'WJets_JESUP',
    'W1Jets_JESUP',
    'W2Jets_JESUP',
    'W3Jets_JESUP',
    'W4Jets_JESUP',
    'TTbar_Madgraph_JESUP',
    'TTbar_Powheg_JESUP',
    'T_t_JESUP',
    'Tbar_t_JESUP',
    'T_tW_JESUP',
    'Tbar_tW_JESUP',
    'T_s_JESUP',
    'Tbar_s_JESUP',
    'WW_JESUP',
    'ZJets_M50_JESUP',
    'ZZ_JESUP',
    #'TT_Mtt_700to1000_JESUP',
    #'TT_Mtt_1000toInf_JESUP',
    'WJets_JESDOWN',
    'W1Jets_JESDOWN',
    'W2Jets_JESDOWN',
    'W3Jets_JESDOWN',
    'W4Jets_JESDOWN',
    'TTbar_Madgraph_JESDOWN',
    'TTbar_Powheg_JESDOWN',
    'T_t_JESDOWN',
    'Tbar_t_JESDOWN',
    'T_tW_JESDOWN',
    'Tbar_tW_JESDOWN',
    'T_s_JESDOWN',
    'Tbar_s_JESDOWN',
    'WW_JESDOWN',
    'ZJets_M50_JESDOWN',
    'ZZ_JESDOWN',    
    #'TT_Mtt_700to1000_JESDOWN',
    #'TT_Mtt_1000toInf_JESDOWN',
    ])    

    if DOJER=='True':
        prefix.extend([
    'WJets_JERUP',
    'W1Jets_JERUP',
    'W2Jets_JERUP',
    'W3Jets_JERUP',
    'W4Jets_JERUP',
    'TTbar_Madgraph_JERUP',
    'TTbar_Powheg_JERUP',
    'T_t_JERUP',
    'Tbar_t_JERUP',
    'T_tW_JERUP',
    'Tbar_tW_JERUP',
    'T_s_JERUP',
    'Tbar_s_JERUP',
    'WW_JERUP',
    'ZJets_M50_JERUP',
    'ZZ_JERUP',
    #'TT_Mtt_700to1000_JERUP',
    #'TT_Mtt_1000toInf_JERUP',
    'WJets_JERDOWN',
    'W1Jets_JERDOWN',
    'W2Jets_JERDOWN',
    'W3Jets_JERDOWN',
    'W4Jets_JERDOWN',
    'TTbar_Madgraph_JERDOWN',
    'TTbar_Powheg_JERDOWN',
    'T_t_JERDOWN',
    'Tbar_t_JERDOWN',
    'T_tW_JERDOWN',
    'Tbar_tW_JERDOWN',
    'T_s_JERDOWN',
    'Tbar_s_JERDOWN',
    'WW_JERDOWN',
    'ZJets_M50_JERDOWN',
    'ZZ_JERDOWN',
    #'TT_Mtt_700to1000_JERDOWN',
    #'TT_Mtt_1000toInf_JERDOWN',
    ])    

    if DOBTAG=='True':
        prefix.extend([
    'WJets_BTAGUP',
    'W1Jets_BTAGUP',
    'W2Jets_BTAGUP',
    'W3Jets_BTAGUP',
    'W4Jets_BTAGUP',
    'TTbar_Madgraph_BTAGUP',
    'TTbar_Powheg_BTAGUP',
    'T_t_BTAGUP',
    'Tbar_t_BTAGUP',
    'T_tW_BTAGUP',
    'Tbar_tW_BTAGUP',
    'T_s_BTAGUP',
    'Tbar_s_BTAGUP',
    'WW_BTAGUP',
    'ZJets_M50_BTAGUP',
    'ZZ_BTAGUP',
    #'TT_Mtt_700to1000_BTAGUP',
    #'TT_Mtt_1000toInf_BTAGUP',
    'WJets_BTAGDOWN',
    'W1Jets_BTAGDOWN',
    'W2Jets_BTAGDOWN',
    'W3Jets_BTAGDOWN',
    'W4Jets_BTAGDOWN',
    'TTbar_Madgraph_BTAGDOWN',
    'TTbar_Powheg_BTAGDOWN',
    'T_t_BTAGDOWN',
    'Tbar_t_BTAGDOWN',
    'T_tW_BTAGDOWN',
    'Tbar_tW_BTAGDOWN',
    'T_s_BTAGDOWN',
    'Tbar_s_BTAGDOWN',
    'WW_BTAGDOWN',
    'ZJets_M50_BTAGDOWN',
    'ZZ_BTAGDOWN',
    #'TT_Mtt_700to1000_BTAGDOWN', 
    #'TT_Mtt_1000toInf_BTAGDOWN',
    ])
    
if DOQCDMC=='True':
    prefix.extend([
    #'QCD_Pt_80to120_Mu',
    #'QCD_Pt_120to170_Mu',
    #'QCD_Pt_170to300_Mu',
    #'QCD_Pt_300to470_Mu',
    #'QCD_Pt_470to600_Mu',
    #'QCD_Pt_600to800_Mu',
    #'QCD_Pt_800to1000_Mu',
    'QCD_Pt_80_170_EM',
    'QCD_Pt_170_250_EM',
    'QCD_Pt_250_350_EM',
    'QCD_Pt_350_EM',
    #'TTbar_MCatNLO'
    ])
    if DOJES == 'True':
        prefix.extend([
        'QCD_Pt_80_170_EM_JESUP',
        'QCD_Pt_170_250_EM_JESUP',
        'QCD_Pt_250_350_EM_JESUP',
        'QCD_Pt_350_EM_JESUP',
        #'TTbar_MCatNLO_JESUP',
        'QCD_Pt_80_170_EM_JESDOWN',
        'QCD_Pt_170_250_EM_JESDOWN',
        'QCD_Pt_250_350_EM_JESDOWN',
        'QCD_Pt_350_EM_JESDOWN',
        #'TTbar_MCatNLO_JESDOWN',
        ])
    if DOJER == 'True':
        prefix.extend([
        'QCD_Pt_80_170_EM_JERUP',
        'QCD_Pt_170_250_EM_JERUP',
        'QCD_Pt_250_350_EM_JERUP',
        'QCD_Pt_350_EM_JERUP',
        #'TTbar_MCatNLO_JERUP',    
        'QCD_Pt_80_170_EM_JERDOWN',
        'QCD_Pt_170_250_EM_JERDOWN',
        'QCD_Pt_250_350_EM_JERDOWN',
        'QCD_Pt_350_EM_JERDOWN',
        #'TTbar_MCatNLO_JERDOWN'    
        ])
    if DOBTAG == 'True':
        prefix.extend([
        'QCD_Pt_80_170_EM_BTAGUP',
        'QCD_Pt_170_250_EM_BTAGUP',
        'QCD_Pt_250_350_EM_BTAGUP',
        'QCD_Pt_350_EM_BTAGUP',
        #'TTbar_MCatNLO_BTAGUP',    
        'QCD_Pt_80_170_EM_BTAGDOWN',
        'QCD_Pt_170_250_EM_BTAGDOWN',
        'QCD_Pt_250_350_EM_BTAGDOWN',
        'QCD_Pt_350_EM_BTAGDOWN',
        #'TTbar_MCatNLO_BTAGDOWN'
        ])

if DOTTBARSYS=='True':
    prefix.extend([
    'TTbar_matchingdown',
    'TTbar_matchingup',
    'TTbar_scaledown',
    'TTbar_scaleup',
    #'TTbar_Powheg_scaleup',
    #'TTbar_Powheg_scaledown',
    ])

###########################################
### Where to save your output root files to
###########################################
dir = []
for i in prefix:
    outdiri = outdir+i
    dir.extend([outdiri])

################################################
### Where is your list of root files to run over
################################################
list = [] 

listnom = [
'Samples_2012/WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_v2_TLBSM_53x_v2_cff.txt',
'Samples_2012/W1JetsToLNu_TuneZ2Star_8TeV_madgraph_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/W2JetsToLNu_TuneZ2Star_8TeV_madgraph_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/W3JetsToLNu_TuneZ2Star_8TeV_madgraph_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/W4JetsToLNu_TuneZ2Star_8TeV_madgraph_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/TT_CT10_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v2_TLBSM_53x_v2_cff.txt',
'Samples_2012/T_t_channel_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/T_s_channel_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/Tbar_s_channel_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/WW_TuneZ2star_8TeV_pythia6_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/ZZ_TuneZ2star_8TeV_pythia6_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/TT_Mtt_700to1000_CT10_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/TT_Mtt_1000toInf_CT10_TuneZ2star_8TeV_powheg_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt' 
    ]

if DONOMINAL=='True':
    list.extend(listnom)
    if DOJES=='True':
        list.extend(listnom)
        list.extend(listnom)
    if DOJER=='True':
        list.extend(listnom)
        list.extend(listnom)
    if DOBTAG=='True':
        list.extend(listnom)
        list.extend(listnom)

listqcd = [
#'Samples_2012/QCD_Pt_80to120_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_120to170_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_170to300_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_300to470_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_470to600_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_600to800_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/QCD_Pt_800to1000_MuEnrichedPt5_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/TT_8TeV_mcatnlo_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
    ] 

if DOQCDMC=='True':
    list.extend(listqcd) #NOMINAL
    if DOJES == 'True':
        list.extend(listqcd) #JESUP
        list.extend(listqcd) #JESDOWN
    if DOJER == 'True':
        list.extend(listqcd) #JERUP
        list.extend(listqcd) #JERDOWN
    if DOBTAG == 'True':
        list.extend(listqcd) #BTAGUP
        list.extend(listqcd) #BTAGDOWN

if DOTTBARSYS=='True':
    list.extend([
'Samples_2012/TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/TT_scaleup_CT10_TuneZ2star_8TeV_powheg_tauola_Summer12_START52_V9_FSIM_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/TT_scaledown_CT10_TuneZ2star_8TeV_powheg_tauola_Summer12_START52_V9_FSIM_v1_TLBSM_53x_v2_cff.txt',
    ])
   

for i in range(len(prefix)):
    print i,' ',prefix[i],' ',dir[i],' ',list[i]

### Write the files you wish to run over for each job    
def get_input(num, list):
    result = '' 
    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    file_count = 0
    for line in file_list:
        if line.find('root')>0:
            file_count=file_count+1
            if file_count>(num-1) and file_count<(num+files_per_job):
                f_name=re.search('.+\'(.+\.root)',line)
                #result=result+'                 \'' + f_name.group(1)+'\',\n'
                result=result+'                 \'dcap:///pnfs/cms/WAX/11' + f_name.group(1)+'\',\n'
    file_list.close()
    #result = result + '                 )\n'
    return result


print str(files_per_job)+' files per job...'

for i in range(len(prefix)):

    j = 1
    nfiles = 1

    ### False by default
    JECUNCERTUP = 'False'
    JECUNCERTDOWN = 'False'
    JERUNCERTUP = 'False'
    JERUNCERTDOWN = 'False'
    BTAGUNCERTUP = 'False'
    BTAGUNCERTDOWN = 'False'

    if (prefix[i].endswith('JESUP')):
        JECUNCERTUP = 'True'
    elif (prefix[i].endswith('JESDOWN')):
        JECUNCERTDOWN = 'True'
    elif (prefix[i].endswith('JERUP')):
        JERUNCERTUP = 'True'
    elif (prefix[i].endswith('JERDOWN')):
        JERUNCERTDOWN = 'True'
    elif (prefix[i].endswith('BTAGUP')):
        BTAGUNCERTUP = 'True'
    elif (prefix[i].endswith('BTAGDOWN')):
        BTAGUNCERTDOWN = 'True'


    ISWJETS = 'False'
    if (prefix[i].startswith('WJets')):
        ISWJETS = 'True'
    if (prefix[i].startswith('W1Jets')):
        ISWJETS = 'True'
    if (prefix[i].startswith('W2Jets')):
        ISWJETS = 'True'
    if (prefix[i].startswith('W3Jets')):
        ISWJETS = 'True'
    if (prefix[i].startswith('W4Jets')):
        ISWJETS = 'True'
                                
    ISTB = 'False'
    if (prefix[i].startswith('T_s') or prefix[i].startswith('Tbar_s')):
        ISTB = 'True'

    ISTT = 'False'
    if (prefix[i].startswith('TT')):
        ISTT = 'True'
        
    print 'CONDOR work dir: '+dir[i]
    #os.system('rm -rf '+dir[i])
    os.system('mkdir -p '+dir[i])

    file_list = open(rel_base+"/src/LJMet/Com/python/"+list[i])
    count = 0
    for line in file_list:
        if line.find('root')>0:
            count = count + 1
    file_list.close()
    #count = count - 1

    print 'File prefix: '+prefix[i]
    print 'Number of input files: '+str(count)

    while ( nfiles <= count ):    

        py_templ_file = open(rel_base+"/src/LJMet/Com/condor/Wprimepython.templ")
        condor_templ_file = open(rel_base+"/src/LJMet/Com/condor/Wprimecondor.templ")
        csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor/Wprimecsh.templ")

        py_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".py","w")
        for line in py_templ_file:
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            line=line.replace('SELECTOR',SELECTOR)
            line=line.replace('INFILES',get_input(nfiles, list[i]))
            line=line.replace('BTAGUNCERTUP',BTAGUNCERTUP)
            line=line.replace('BTAGUNCERTDOWN',BTAGUNCERTDOWN)
            line=line.replace('JECUNCERTUP',JECUNCERTUP)
            line=line.replace('JECUNCERTDOWN',JECUNCERTDOWN)
            line=line.replace('JERUNCERTUP',JERUNCERTUP)
            line=line.replace('JERUNCERTDOWN',JERUNCERTDOWN)
            line=line.replace('EVENTSTOPROCESS',str(-1))
            line=line.replace('ISMCSAMPLE',ISMCSAMPLE)
            line=line.replace('JSONFILE',MYJSON)
            line=line.replace('DOLASERCALFILT',DOLASERCALFILT)
            line=line.replace('MCTRIGGEREL',MCTRIGGEREL) 
            line=line.replace('MCTRIGGERMU',MCTRIGGERMU) 
            line=line.replace('DATATRIGGEREL',DATATRIGGEREL) 
            line=line.replace('DATATRIGGERMU',DATATRIGGERMU)
            line=line.replace('MINTIGHTMUON',MINTIGHTMUON)
            line=line.replace('MINTIGHTELECTRON',MINTIGHTELECTRON)
            line=line.replace('MINTIGHTLEPTON',MINTIGHTLEPTON)
            line=line.replace('MAXTIGHTLEPTON',MAXTIGHTLEPTON)
            line=line.replace('LEADINGJETPT',LEADINGJETPT)
            line=line.replace('ISWJETS',ISWJETS)
            line=line.replace('ISTB',ISTB)
            line=line.replace('ISTT',ISTT)
            line=line.replace('USEBESTTOP',USEBESTTOP)
            line=line.replace('DO53XJEC',DO53XJEC)
            py_file.write(line)
        py_file.close()

        condor_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".condor","w")
        for line in condor_templ_file:
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            condor_file.write(line)
        condor_file.close()

        csh_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".csh","w")
        for line in csh_templ_file:
            line=line.replace('CMSSWBASE',rel_base)
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            line=line.replace('FWLITEANALYZER',FWLiteAnalyzer)
            csh_file.write(line)
        csh_file.close()

        os.system('chmod u+x '+dir[i]+'/'+prefix[i]+'_'+str(j)+'.csh')
        os.system('cd '+dir[i]+'; condor_submit '+prefix[i]+'_'+str(j)+'.condor; cd -')
        j = j + 1
        nfiles = nfiles + files_per_job
        py_templ_file.close()
        condor_templ_file.close()
        csh_templ_file.close()

    print  str(j-1)+' jobs submitted'
    
