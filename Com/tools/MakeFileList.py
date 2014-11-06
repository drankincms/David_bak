#!/usr/bin/python

import os
import re
import fileinput
import glob



dir = [   '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_800_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_800_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_800_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_900_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_900_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1000_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1000_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1000_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1100_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1100_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1100_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1200_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1200_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1200_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1300_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1300_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1300_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1400_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1400_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1400_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1500_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1500_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1500_LeftWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1600_RightWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1600_MixRLWprime',
          '/mnt/hadoop/data/se/store/user/msegala/Wprime/WprimeToTB_M_1600_LeftWprime'
          ]


prefix = ['WprimeToTB_M_800_RightWprime',
          'WprimeToTB_M_800_MixRLWprime',
          'WprimeToTB_M_800_LeftWprime',
          'WprimeToTB_M_900_RightWprime',
          'WprimeToTB_M_900_MixRLWprime',
          'WprimeToTB_M_1000_RightWprime',
          'WprimeToTB_M_1000_MixRLWprime',
          'WprimeToTB_M_1000_LeftWprime',
          'WprimeToTB_M_1100_RightWprime',
          'WprimeToTB_M_1100_MixRLWprime',
          'WprimeToTB_M_1100_LeftWprime',
          'WprimeToTB_M_1200_RightWprime',
          'WprimeToTB_M_1200_MixRLWprime',
          'WprimeToTB_M_1200_LeftWprime',
          'WprimeToTB_M_1300_RightWprime',
          'WprimeToTB_M_1300_MixRLWprime',
          'WprimeToTB_M_1300_LeftWprime',
          'WprimeToTB_M_1400_RightWprime',
          'WprimeToTB_M_1400_MixRLWprime',
          'WprimeToTB_M_1400_LeftWprime',
          'WprimeToTB_M_1500_RightWprime',
          'WprimeToTB_M_1500_MixRLWprime',
          'WprimeToTB_M_1500_LeftWprime',
          'WprimeToTB_M_1600_RightWprime',
          'WprimeToTB_M_1600_MixRLWprime',
          'WprimeToTB_M_1600_LeftWprime'
          ]



output = '/data/users/segala/test/CMSSW_4_2_4_patch2/src/LJMet/Com/python/Pat_Layer_1_2011_Summer_MC_Wprime/'

j = 0 
for i in dir:

    outputFile = output + prefix[j] + '.py'

    f = open(outputFile,'w')
    
    for infile in glob.glob( os.path.join(i, '*.root') ):
        f.write(infile + '\n')
        
        
    f.close()
    j = j +1
