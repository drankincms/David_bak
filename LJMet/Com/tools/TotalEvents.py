import os
import re
import glob
import sys

#dir = '/uscms_data/d2/mike1886/Top/Wprime_withWeights_Nov15/WJets_scaledown_2/'

#dir = '/uscms_data/d2/lpcbtag/segala/TOP/WprimeFastSim_Dec22/WprimeToTB_M_1600_LeftWprime'
#dir = '/uscms_data/d2/samirg28/PUC/Wprime/Muon_Files_NewTest/Wei_Data_Run2011B_958pb/'
#Wei_Data_Run2011A_2418pb/  Wei_Data_Run2011B_1796pb/  Wei_Data_Run2011B_958pb/

dir = '/uscms_data/d2/mike1886/LJMet/ana_output/SingleMu_2012A_598pb'

print 'For file: ',dir  

### Open your input file and get the cut flow

total = 0
i = 1
for infile in glob.glob( os.path.join(dir, '*.stdout') ):
    #print "current file is: " + infile
    print 'i = ',i
    file = open(infile)
    for line in file:

        if line.find('No selection')>0:
            f_name=re.search('(?<=selection    )\w+',line)
            if f_name:
                print f_name.group(0)
                total = total + int(f_name.group(0))
            else:
                f_name=re.search('(?<=selection   )\w+',line)
                if f_name:
                    print f_name.group(0)
                    total = total + int(f_name.group(0))
                else:
                    f_name=re.search('(?<=selection  )\w+',line)
                    if f_name:
                        print f_name.group(0)
                        total = total + int(f_name.group(0))
                    else:
                        f_name=re.search('(?<=selection     )\w+',line)
                        if f_name:
                            print f_name.group(0)
                            total = total + int(f_name.group(0))
                        else:
                            f_name=re.search('(?<=selection      )\w+',line)
                            if f_name:
                                print f_name.group(0)
                                total = total + int(f_name.group(0))
                            else:
                                f_name=re.search('(?<=selection       )\w+',line)
                                if f_name:
                                    print f_name.group(0)
                                    total = total + int(f_name.group(0))
    i = i + 1
        

print 'Total = ',total


#     0 :         No selection    1667889
