#!/usr/bin/python

import os
import re
import fileinput


files_per_job =7

rel_base = os.environ['CMSSW_BASE']


#################################################
### What names to give to your output root files
#################################################

#prefix = 'WJets'
#prefix = 'WJets_Pythia'
prefix = 'WJets_PileUp'
#prefix = 'WJets_Madgraph'
#prefix = 'Inclusive_QCD'
#prefix = 'TTbar'
#prefix = 'Wtaunu'
#prefix = 'Ztautau'
#prefix = 'Zmumu'

#prefix = 'DataSample'
#prefix = 'MinBias_pythia8'


###########################################
### Where to save your output root files to
###########################################

#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ReverseIso/WJets'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ReverseIso/QCD'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ReverseIso/Data_3pb'




####### FOR JINST

#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/Full_DataSample_8pb'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/TTbar'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/Zmumu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/Wtaunu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/W_Jets'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/Inclusive_QCD'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_10252010/Ztautau'

## For 38 samples
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/Full_DataSample_35pb'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/W_Jets_WithHLT'
dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/W_Jets_WithPileUp'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/Zmumu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/Zmumu_WithPileUp'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/TTBar'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/WTauNu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/WTauNu_WithPileUp'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/ZTauTau'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/ZTauTau_WithPileUp'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/QCD'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForJINST_38Samples_1122010/QCD_WithPileUp'




#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_812010/Full_DataSample_3pb'

#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_812010/MinBias_3GeV'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/TTbar'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/Zmumu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/Wtaunu'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/W_Jets'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/Inclusive_QCD'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_9132010/Ztautau'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_812010/W_Pythia_GenW'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_712010/W_PileUp_PVZvalues'
#dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_812010/W_Madgraph'

################################################
### Where is your list of root files to run over
################################################


MultpleFiles = 0

if MultpleFiles == 1:

   ####list = 'Pat_Layer1_June14th/MinBias_712010.py'
   ####list1 = 'Pat_Layer1_June14th/Run2010A_712010.py'
   ####list2 = 'Pat_Layer1_June14th/Run2010A_PromptReco_7152010.py'
    
#### the top 2 are for 246 nb-1
    #list = 'Pat_Layer1_June14th/7212010_July16thReReco.py'
    #list1 = 'Pat_Layer1_June14th/Run2010A_PromptReco_7212010.py'
    #####list1 = 'Pat_Layer1_June14th/Run2010A_PromptReco_816.py'


#### For 8 pb-1 sample
    #list = 'Pat_Layer1_Mu_Run2010A_Sep17ReReco_v2/Mu_Run2010A_Sep17ReReco_v2_10192010.py'
    #list1 = 'Pat_Layer1_Mu_Run2010A_Sep17ReReco_v2/Mu_Run2010B_PromptReco_v2_20Oct2010.py'

    
#### For 35 pb-1 sample
    list = 'Pat_Layer1_Mu_Run2010A_Sep17ReReco_v2/Mu_Run2010A_Sep17ReReco_v2_10192010.py'
    list1 = 'Pat_Layer1_MC_38_Samples/Mu_Run2010B_PromptReco_v2_10Nov2010.py'
    
    






    #list = 'Pat_Layer1_MC_722010/WplusToMuNu_CTEQ66_powheg_Summer10_START36_V9_S09_v1_722010.py'
    #list1 = 'Pat_Layer1_MC_722010/WminusToMuNu_CTEQ66_powheg_Summer10_START36_V9_S09_v1_722010.py'

    #list = 'Pat_Layer1_MC_722010/WplusToTauNu_CTEQ66_powheg_Summer10_START36_V9_S09_v2_722010.py'
    #list1 = 'Pat_Layer1_MC_722010/WminusToTauNu_CTEQ66_powheg_Summer10_START36_V9_S09_v2_722010.py'

    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    file_list1= open(rel_base+"/src/LJMet/Com/python/"+list1)
    #file_list2= open(rel_base+"/src/LJMet/Com/python/"+list2)
       

    ALLfiles = []
    ALLfiles.append(file_list)
    ALLfiles.append(file_list1)
    #ALLfiles.append(file_list2)


    count = 0
    newfile = open("TotalListFile.py","w")
    for i in range(2):
        for line in ALLfiles[i]:
            if line.find('root')>0:
                newfile.write(line)
                count = count + 1
    newfile.close()
    print "Total Entries = ",count
    count = count - 1

else:  
    #list = 'Pat_Layer1_MC_722010/TTbar_Summer10_START36_V9_S09_v1_722010.py'
    #list = 'Pat_Layer1_MC_722010/Zmumu_M20_CTEQ66_powheg_Summer10_START36_V9_S09_v2_722010.py'
    #list = 'Pat_Layer1_MC_722010/InclusiveMu15_Summer10_START36_V9_S09_v1_722010.py'
    #list = 'Pat_Layer1_MC_722010/Ztautau_M20_CTEQ66_powheg_Summer10_START36_V9_S09_v1_742010.py'
    #list = 'Pat_Layer1_MC_722010/Wmunu_Summer10_START36_V9_S09_v2_7122010.py'
    #list = 'PileUP_Pat_Layer1/W_Jets_NoPileUp_782010.py'
    #list = 'PileUP_Pat_Layer1/W_Jets_WithPileUp_782010.py'
    
    #For 38 samples
    #list = 'Pat_Layer1_MC_38_Samples/WToMuNu_TuneZ2_7TeV_pythia6_Fall10_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/WToMuNu_TuneZ2_7TeV_pythia6_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/DYToMuMu_M_20_CT10_TuneZ2_7TeV_poweg_pythia6_Fall10_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/DYToMuMu_M_20_CT10_TuneZ2_7TeV_poweg_pythia6_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'    
    #list = 'Pat_Layer1_MC_38_Samples/TT_TuneZ2_7TeV_pythia6_tauola_Fall10_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/DYToTauTau_M_20_TuneZ2_7TeV_pythia6_tauola_Fall10_START38_V12_v1.py' 
    #list = 'Pat_Layer1_MC_38_Samples/DYToTauTau_M_20_TuneZ2_7TeV_pythia6_tauola_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/WToTauNu_TuneZ2_7TeV_pythia6_tauola_Fall10_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/WToTauNu_TuneZ2_7TeV_pythia6_tauola_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6_Fall10_START38_V12_v1.py'
    #list = 'Pat_Layer1_MC_38_Samples/QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'


    #For 38 samples with corrected pf type1

    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)

    list = 'Pat_Layer1_MC_38_Samples_FixedPfMET/DYToMuMu_M_20_CT10_TuneZ2_7TeV_poweg_pythia6_Fall10_E7TeV_ProbDist_2010Data_BX156_START38_V12_v1.py'

    count = 0
    for line in file_list:
        if line.find('root')>0:
            count = count + 1
    file_list.close()
    count = count - 1




print 'CONDOR work dir: '+dir
print 'File prefix: '+prefix
print 'Number of input files: '+str(count)
print str(files_per_job)+' files per job...'


os.system('rm -rf '+dir)
os.system('mkdir -p '+dir)


### Write the files you wish to run over for each job    
def get_input(num):
    result =          'process.inputs = cms.PSet (\n'
    result = result + '    fileNames = cms.vstring(\n'
    if MultpleFiles == 1:
        file_list = open(rel_base+"/src/LJMet/Com/condor/TotalListFile.py")
    else:    
        file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    file_count = 0
    for line in file_list:
        if line.find('root')>0:
            file_count=file_count+1
            if file_count>(num-1) and file_count<(num+files_per_job):
                f_name=re.search('.+\'(.+\.root)',line)
                result=result+'    \'' + f_name.group(1)+'\',\n'
    file_list.close()
    result = result + '    )\n'
    result = result + ')\n'     
    return result

nfiles = 1
i = 1


### Which MET collection do you want to run over
met = 'patMETs'

### What is the name of your FWLite Analyzer
FWLiteAnalyzer = 'PatFwliteALLMetAnalyzer_Selector'

# Create and run all python files

while ( nfiles <= count ):    
   
    py_templ_file = open(rel_base+"/src/LJMet/Com/condor/python.templ")
    condor_templ_file = open(rel_base+"/src/LJMet/Com/condor/condor.templ")
    csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor/csh.templ")
    
    py_file = open(dir+"/"+prefix+"_"+str(i)+".py","w")
    for line in py_templ_file:
        line=line.replace('DIRECTORY',dir)
        line=line.replace('PREFIX',prefix)
        line=line.replace('JOBID',str(i))
        line=line.replace('CFISOURCE',get_input(nfiles))
        line=line.replace('EVENTSTOPROCESS',str(-1))
        line=line.replace('METSOURCE',met)
        py_file.write(line)
    py_file.close()
    
    condor_file = open(dir+"/"+prefix+"_"+str(i)+".condor","w")
    for line in condor_templ_file:
        line=line.replace('DIRECTORY',dir)
        line=line.replace('PREFIX',prefix)
        line=line.replace('JOBID',str(i))
        condor_file.write(line)
    condor_file.close()

    csh_file = open(dir+"/"+prefix+"_"+str(i)+".csh","w")
    for line in csh_templ_file:
        line=line.replace('CMSSWBASE',rel_base)
        line=line.replace('DIRECTORY',dir)
        line=line.replace('PREFIX',prefix)
        line=line.replace('JOBID',str(i))
        line=line.replace('FWLITEANALYZER',FWLiteAnalyzer)
        csh_file.write(line)
    csh_file.close()

    os.system('chmod u+x '+dir+'/'+prefix+'_'+str(i)+'.csh')
    os.system('condor_submit '+dir+'/'+prefix+'_'+str(i)+'.condor')
    i = i + 1
    nfiles = nfiles + files_per_job
    py_templ_file.close()
    condor_templ_file.close()
    csh_templ_file.close()
    
print  str(i-1)+' jobs submitted'

    
 
