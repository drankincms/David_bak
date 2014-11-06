#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 50

rel_base = os.environ['CMSSW_BASE']

###########################################
### Where to save your output root files to
###########################################
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

### Is this a TB sample?
ISTB = 'True'
ISTT = 'False'

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
DONOM = 'True'
DOJES = 'True'
DOJER = 'True'
DOBTAG = 'True'


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

Rightmasses = ['Wprime800Right','Wprime900Right','Wprime1000Right','Wprime1100Right','Wprime1200Right','Wprime1300Right','Wprime1400Right','Wprime1500Right','Wprime1600Right','Wprime1700Right','Wprime1800Right','Wprime1900Right','Wprime2000Right','Wprime2100Right','Wprime2200Right','Wprime2300Right','Wprime2400Right','Wprime2500Right','Wprime2600Right','Wprime2700Right','Wprime2800Right','Wprime2900Right','Wprime3000Right']

Rightlist = [ 
'Samples_2012/SingletopWprime_M_800_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_900_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1000_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1100_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1200_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1300_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1400_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1500_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1600_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1700_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1800_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1900_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2000_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2100_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2200_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2300_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2400_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2500_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2600_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2700_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2800_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2900_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_3000_right_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
]

Leftmasses = ['Wprime800Left','Wprime900Left','Wprime1000Left','Wprime1100Left','Wprime1200Left','Wprime1300Left','Wprime1400Left','Wprime1500Left','Wprime1600Left','Wprime1700Left','Wprime1800Left','Wprime1900Left','Wprime2000Left','Wprime2100Left','Wprime2200Left','Wprime2300Left','Wprime2400Left','Wprime2500Left','Wprime2600Left','Wprime2700Left','Wprime2800Left','Wprime2900Left','Wprime3000Left']

Leftlist = [ 
'Samples_2012/SingletopWprime_M_800_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_900_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1000_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1100_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1200_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1300_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1400_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1500_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1600_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1700_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1800_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1900_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2000_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2100_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2200_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2300_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2400_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2500_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2600_left_TuneZ2star_8TeV_comphep_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2700_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2800_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2900_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_3000_left_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
]


Mixmasses = ['Wprime800Mix','Wprime900Mix','Wprime1000Mix','Wprime1100Mix','Wprime1200Mix','Wprime1300Mix','Wprime1400Mix','Wprime1500Mix','Wprime1600Mix','Wprime1700Mix','Wprime1800Mix','Wprime1900Mix','Wprime2000Mix','Wprime2100Mix','Wprime2200Mix','Wprime2300Mix','Wprime2400Mix','Wprime2500Mix','Wprime2600Mix','Wprime2700Mix','Wprime2800Mix','Wprime2900Mix']

Mixlist = [ 
'Samples_2012/SingletopWprime_M_800_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_900_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1000_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1100_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1200_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1300_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1400_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1500_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1600_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1700_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1800_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_1900_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2000_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2100_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2200_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2300_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2400_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2500_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2600_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2700_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2800_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingletopWprime_M_2900_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
#'Samples_2012/SingletopWprime_M_3000_mixed_TuneZ2star_8TeV_comphep_StoreResults_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_cff.txt',
]

prefix = []
list = []

i=0

print Rightlist[i]
for rm in Rightmasses:
    if DONOM == 'True':
        prefix.extend([rm])
        list.extend([Rightlist[i]])
    if DOJES == 'True':
        prefix.extend([rm+'_JESUP'])
        prefix.extend([rm+'_JESDOWN'])
        list.extend([Rightlist[i]])
        list.extend([Rightlist[i]])
    if DOJER == 'True':
        prefix.extend([rm+'_JERUP'])
        prefix.extend([rm+'_JERDOWN'])
        list.extend([Rightlist[i]])
        list.extend([Rightlist[i]])
    if DOBTAG == 'True':
        prefix.extend([rm+'_BTAGUP'])
        prefix.extend([rm+'_BTAGDOWN'])
        list.extend([Rightlist[i]])
        list.extend([Rightlist[i]])
    i+=1 

i=0
for lm in Leftmasses: 
    if DONOM == 'True':
        prefix.extend([lm])
        list.extend([Leftlist[i]])
    if DOJES == 'True':
        prefix.extend([lm+'_JESUP'])
        prefix.extend([lm+'_JESDOWN'])
        list.extend([Leftlist[i]])
        list.extend([Leftlist[i]])
    if DOJER == 'True':
        prefix.extend([lm+'_JERUP'])
        prefix.extend([lm+'_JERDOWN'])
        list.extend([Leftlist[i]])
        list.extend([Leftlist[i]])
    if DOBTAG == 'True':
        prefix.extend([lm+'_BTAGUP'])
        prefix.extend([lm+'_BTAGDOWN'])
        list.extend([Leftlist[i]])
        list.extend([Leftlist[i]])
    i+=1
i=0
for mm in Mixmasses:
    if DONOM == 'True':
        prefix.extend([mm])
        list.extend([Mixlist[i]])
    if DOJES == 'True':
        prefix.extend([mm+'_JESUP'])
        prefix.extend([mm+'_JESDOWN'])
        list.extend([Mixlist[i]])
        list.extend([Mixlist[i]])
    if DOJER == 'True':
        prefix.extend([mm+'_JERUP'])
        prefix.extend([mm+'_JERDOWN'])
        list.extend([Mixlist[i]])
        list.extend([Mixlist[i]])
    if DOBTAG == 'True':
        prefix.extend([mm+'_BTAGUP'])
        prefix.extend([mm+'_BTAGDOWN'])
        list.extend([Mixlist[i]])
        list.extend([Mixlist[i]])
    i+=1
 
dir = []
for i in prefix:
    dir.extend([outdir+i])

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


    print i,' CONDOR work dir: '+dir[i],' list ',list[i]
    os.system('rm -rf '+dir[i])
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
    
