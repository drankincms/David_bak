import os
import re
import glob
import sys

#dir = '/uscms_data/d2/mike1886/Top/Wprime_withWeights_Dec7/'
#dir = '/uscms_data/d2/mike1886/Top/Wprime_withWeights_Dec18/'


#dir = '/uscms_data/d2/lpcbtag/segala/TOP/WprimeOneJetRequirment_Dec22/'
#dir = '/uscms_data/d2/lpcbtag/segala/TOP/Wprime_Dec22/'

dir = '/uscms_data/d2/mike1886/LJMet/ana_output/'



'''
list = ['Data_Run2011A_2418pb',
        'Data_Run2011B_1796pb',
        'Data_Run2011B_958pb'
        ]

prefix = [ 'DataSample', 
           'DataSample',
           'DataSample'
           ]
'''



'''
list = ['TTbar',
        'WJets',
        'WW',
        'ZJets',
        'QCD',
        'T_t',
        'Tbar_t',
        'T_tW',
        'Tbar_tW',
        'T_s',
        'Tbar_s'
        ]


prefix = [ 'TTbar', 
           'WJets', 
           'WW', 
          'ZJets', 
          'Inclusive_QCD', 
          'T_t', 'Tbar_t','T_tW', 'Tbar_tW','T_s', 'Tbar_s'
          ]
'''



list = ['TTbar_JER_UP',
        'WJets_JER_UP',
        #'WW_JER_UP',
        #'ZJets_JER_UP',
        #'QCD_JER_UP',
        #'T_t_JER_UP',
        'Tbar_t_JER_UP',
        'T_tW_JER_UP',
        'Tbar_tW_JER_UP',
        #'T_s_JER_UP',
        #'Tbar_s_JER_UP',
        'TTbar_JER_DOWN',
        'WJets_JER_DOWN',
        #'WW_JER_DOWN',
        #'ZJets_JER_DOWN',
        #'QCD_JER_DOWN',
        #'T_t_JER_DOWN',
        'Tbar_t_JER_DOWN',
        'T_tW_JER_DOWN',
        'Tbar_tW_JER_DOWN'
        #'T_s_JER_DOWN',
        #'Tbar_s_JER_DOWN'
        ]

prefix = [ 'TTbar', 
           'WJets', 
           #'WW', 
           #'ZJets', 
           #'Inclusive_QCD', 
           #'T_t', 
           'Tbar_t','T_tW', 'Tbar_tW',
           #'T_s', 'Tbar_s',
           'TTbar', 
           'WJets', 
           #'WW', 
           #'ZJets', 
           #'Inclusive_QCD', 
           #'T_t', 
           'Tbar_t','T_tW', 'Tbar_tW'
           #'T_s', 'Tbar_s'
           ]









'''

list = ['TTbar_matchingup',
        'TTbar_matchingdown',
        'TTbar_scaleup',
        'TTbar_scaledown',
        'WJets_matchingup',
        'WJets_matchingdown',
        'WJets_scaleup',
        'WJets_scaledown'
        ]

prefix = ['TTbar_matchingup',
          'TTbar_matchingdown',
          'TTbar_scaleup',
          'TTbar_scaledown',
          'WJets_matchingup',
          'WJets_matchingdown',
          'WJets_scaleup',
          'WJets_scaledown'
          ]
'''




j = 0
for i in list:
    direc = dir + i

    os.chdir( direc )
    os.system('hadd ' + prefix[j] + '.root ' + prefix[j] +'_*' +'root')
    #os.system('rm All.root')


    j += 1
