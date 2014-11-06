###########################################################################################
########
######## This script locates all the condor files produced and adds up thier cut flow table
########
###########################################################################################

import os
import re


dir = '/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/MinimumBias_SD_MU/CaloLooseQaulityCuts_WithTrees_JetBins/'
file = 'MinimumBias_Commissioning10_SD_MU_'
suffix = '.stdout'
i = 1

file_list = open(dir + file + str(i) + suffix)

### Open your input file and get the cut flow


for line in file_list:
    if line.find('Trigger cuts')>0:
        print 'YES I FOUND IT'
        f_name=re.search('(Trigger cuts)',line)
        print f_name.group(0)
file_list.close()



#Trigger cuts
