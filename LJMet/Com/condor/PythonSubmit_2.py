#!/usr/bin/python

import os
import re
import fileinput


files_per_job = 20

rel_base = os.environ['CMSSW_BASE']


#################################################
### What names to give to your output root files
#################################################

#prefix = 'MinimumBias_Commissioning10_May6thPDSkim_GOODCOLL_v1'
#prefix = 'MinimumBias_Commissioning10_SD_MU'
#prefix = 'MinBias_pythia8_Spring10-START3X_V1_5202010'

prefix = 'Full_DataSample'
#prefix = 'Wtaunu'
#prefix = 'Inclusive_QCD'


###########################################
### Where to save your output root files to
###########################################


dir = '/uscms_data/d2/kukarzev/tmp/EWK_MET/MetAnalZ/Data/'
#dir = '/uscms_data/d2/kukarzev/tmp/EwkMETCommissioning/FWLiteAnalyzerOuput/Wtaunu_Spring10/'
#dir = '/uscms_data/d2/kukarzev/tmp/EwkMETCommissioning/FWLiteAnalyzerOuput/Inclusive_QCD/'


################################################
### Where is your list of root files to run over
################################################


MultpleFiles = 1

if MultpleFiles == 1:

    list = 'Pat_Layer1_MinimumBias_Commissioning10_SD_MU/Pat_Layer1_MinimumBias_Commissioning10_SD_MU_672010_May6thPDSkim.py'
    list1 = 'Pat_Layer1_MinimumBias_Commissioning10_SD_MU/Pat_Layer1_MinimumBias_Commissioning10_SD_MU_672010.py'
    list2 = 'Pat_Layer1_Run2010A_MD_SD/Pat_Layer1_Run2010A_MD_SD_672010.py'


    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    file_list1= open(rel_base+"/src/LJMet/Com/python/"+list1)
    file_list2= open(rel_base+"/src/LJMet/Com/python/"+list2)


    #list = 'Pat_Layer1_MC_Samples/Wmunu_Wminus_powheg_Spring10_START3X_V26_S09_v1.py'
    #list1 = 'Pat_Layer1_MC_Samples/Wmunu_Wplus_powheg_Spring10_START3X_V26_S09_v1.py'
    

    #file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    #file_list1= open(rel_base+"/src/LJMet/Com/python/"+list1)
       

    ALLfiles = []
    ALLfiles.append(file_list)
    ALLfiles.append(file_list1)
    ALLfiles.append(file_list2)


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
    list = 'Pat_Layer1_MC_Samples/InclusiveMu15_Spring10_START3X_V26_S09_v1.py'
    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)

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
FWLiteAnalyzer = 'MetAnalZ_PatFwlite'

# Create and run all python files

while ( nfiles < count ):    
   
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

    
 
