#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 30

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
ISMCSAMPLE = 'False'

###Is WJets?
ISWJETS = 'False'

###53x JEC?
DO53XJEC = 'True'

###Is TB?
ISTB = 'False'
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

### JSON file to use
MYJSON = "''"

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

prefix = [
'SingleMu_13Jul2012A',
'SingleMu_06Aug2012A',
'SingleMu_13Jul2012B',
'SingleMu_24Aug2012C',
'SingleMu_Prompt2012C1',
'SingleMu_Prompt2012C2',
'SingleMu_Prompt2012D1',
'SingleMu_Prompt2012D2',
'SingleMu_11Dec2012C',
'SingleElectron_13Jul2012A',
'SingleElectron_06Aug2012A',
'SingleElectron_13Jul2012B',
'SingleElectron_24Aug2012C',
'SingleElectron_Prompt2012C1',
'SingleElectron_Prompt2012C2',
'SingleElectron_Prompt2012D1',
'SingleElectron_Prompt2012D2',
'SingleElectron_11Dec2012C',
]

dir = [
outdir+'SingleMu_13Jul2012A',
outdir+'SingleMu_06Aug2012A',
outdir+'SingleMu_13Jul2012B',
outdir+'SingleMu_24Aug2012C',
outdir+'SingleMu_Prompt2012C1',
outdir+'SingleMu_Prompt2012C2',
outdir+'SingleMu_Prompt2012D1',
outdir+'SingleMu_Prompt2012D2',
outdir+'SingleMu_11Dec2012C',
outdir+'SingleElectron_13Jul2012A',
outdir+'SingleElectron_06Aug2012A',
outdir+'SingleElectron_13Jul2012B',
outdir+'SingleElectron_24Aug2012C',
outdir+'SingleElectron_Prompt2012C1',
outdir+'SingleElectron_Prompt2012C2',
outdir+'SingleElectron_Prompt2012D1',
outdir+'SingleElectron_Prompt2012D2',
outdir+'SingleElectron_11Dec2012C',
]


################################################
### Where is your list of root files to run over
################################################
list = [ 
'Samples_2012/SingleMu_StoreResults_Run2012A_13Jul2012_v1_TLBSM_53x_v2_jsonfix_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012A_recover_06Aug2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012B_13Jul2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012C_24Aug2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012C_PromptReco_v2_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012C_PromptReco_v2_TLBSM_53x_v2_extension_v1_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012D_PromptReco_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleMu_StoreResults_Run2012D_PromptReco_v1_TLBSM_53x_v2_extension_v2_cff.txt',
'Samples_2012/SingleMu_Run2012C_EcalRecover_11Dec2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012A_13Jul2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012A_recover_06Aug2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012B_13Jul2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012C_24Aug2012_v1_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012C_PromptReco_v2_TLBSM_53x_v2_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012C_PromptReco_v2_TLBSM_53x_v2_extension_v1_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012D_PromptReco_v1_TLBSM_53x_v2_bugfix_cff.txt',
'Samples_2012/SingleElectron_StoreResults_Run2012D_PromptReco_v1_TLBSM_53x_v2_extension_v1_cff.txt',
'Samples_2012/SingleElectron_Run2012C_EcalRecover_11Dec2012_v1_TLBSM_53x_v2_cff.txt',
] 

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

    MINTIGHTMUON = '0'
    MINTIGHTELECTRON ='0'

    if (prefix[i].startswith('SingleMu')):
        MINTIGHTMUON = '1'
    elif (prefix[i].startswith('SingleElectron')):
        MINTIGHTELECTRON = '1'

    MYJSON = "''"
    if (prefix[i].find('13Jul2012')>0):
        MYJSON = "'../data/json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt'"
        DOLASERCALFILT='True'
    if (prefix[i].find('06Aug2012')>0):
        MYJSON = "'../data/json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt'"
        DOLASERCALFILT='False'
    if (prefix[i].find('24Aug2012')>0):
        MYJSON = "'../data/json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'"
        DOLASERCALFILT='False'
    if (prefix[i].find('Prompt2012C')>0):
        MYJSON = "'../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_198941-203709.txt'"
        DOLASERCALFILT='False'
    if (prefix[i].find('Prompt2012D')>0):
        MYJSON = "'../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_203894-208686.txt'"
        DOLASERCALFILT='False'
    if (prefix[i].find('11Dec2012')>0):
        MYJSON = "'../data/json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt'"
        DOLASERCALFILT='False'

        
    print 'CONDOR work dir: '+dir[i]
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
    print 'JSON: '+MYJSON

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
        #os.system('cd '+dir[i]+'; condor_submit '+prefix[i]+'_'+str(j)+'.condor; cd -')
                
        j = j + 1
        nfiles = nfiles + files_per_job
        py_templ_file.close()
        condor_templ_file.close()
        csh_templ_file.close()

    print  str(j-1)+' jobs submitted'

    
