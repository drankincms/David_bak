##############################################################################################
##############################################################################################
######
######  This script looks through all the stdout files in a given directory and finds the 
######  original number of events processed in crab
######
######  To use: python CrabEvents.py /path/to/directory
######
######
##############################################################################################
##############################################################################################

import os
import re
import glob
import sys


dir = sys.argv[1] + '/res/'
print 'For file: ',dir  


runs = [1,7,9,10,12,14,18,21,25,26,32,35,37,41,43,44,46,47,54,55,59,60,63,67,68,69,74,75,77,78,79,80,81,83,85,88,93,94,97,106,107,109,110,113,114,119,130,134,136,137,138,140,145,149,154,155,158,159,160,166,172,174] 

total = 0
for infile in glob.glob( os.path.join(dir, '*.stdout') ):
    print "current file is: " + infile

#for i in runs:
    #infile = dir + 'CMSSW_' + str(i) + '.stdout'
    #print "current file is: " + infile
    file_list = open(infile)
    for line in file_list:
        if line.find('total =')>0:
            f_name=re.search('(?<=TrigReport Events total = )\w+',line)
            print f_name.group(0)
            total = total + int(f_name.group(0))


print 'Total events = ',total

