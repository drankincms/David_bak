#####################################################################################################################################################################################
######
###### This python script looks at the output directory created by condor. All files are searched for thier exit codes and if any errors have been found those files are resubmitted
######
###### To run: python CondorResubmitScript.py
######
###### Created by: Michael Segala 
######
#####################################################################################################################################################################################


import os
import re


Dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_712010/Full_DataSample_15GeV' 

files = 0
rootfiles = 0
goodfiles = 0
badfiles = 0

badfilelist = []
badfileerror = []
badrootfiles = []

for f in os.listdir(Dir):
    if f.find('.root')>0:
        rootfiles = rootfiles + 1
        #print f 
    if f.find('condor.log')>0:
        files = files + 1
        file = open(Dir+"/"+f)
        for line in file:
            if line.find('return value')>0:
                f_name=re.search('(return value 0)',line)
                if f_name:
                    goodfiles = goodfiles + 1
                                
                f_name_bad=re.search('return value [^0]',line)
                if f_name_bad:
                    badfiles = badfiles + 1
                    badfilelist.append(f)
                    badfileerror.append(f_name_bad.group(0))
        file.close()


unseccues = files - rootfiles

### This is needed if you have combined all your root files already
if unseccues == -1:
    unseccues = 0

print 'Total number of files ran over = ', files
print 'Number of files not returning a root file = ', unseccues
print 'Number of files exited successfully = ', goodfiles
print 'Number of files exited with an error = ', badfiles
print '-------------------------------------'

for i in range(len(badfilelist)):
    print 'File: ' + badfilelist[i] + ' exited with ' + badfileerror[i]
