import os
import re
import glob
import sys

indir = '/uscms_data/d2/dsperka/8TeV/Samples/29May/'
outdir = '/uscms_data/d2/dsperka/8TeV/Samples/29May_All/'

masses = ['800','900','1000','1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200','2300','2400','2500','2600','2700','2800','2900','3000']
couplings = ['Left','Mix']
#couplings = ['Right']


#if not os.path.exists(outdir):
#    os.mkdir(outdir)

for mass in masses:
    for coup in couplings:  
 
        indirec = indir+'Wprime'+mass+coup
        #print 'hadd ' + outdir + 'Wprime'+mass+coup+'.root' + indirec +'/*.root'
        os.system('hadd '+outdir+'Wprime'+mass+coup+'.root '+indirec+'/*.root')

    
systematics = ['JESUP','JERUP','BTAGUP','JESDOWN','JERDOWN','BTAGDOWN']

for mass in masses:
    for coup in couplings:  
        for sys in systematics:

            indirec = indir+'Wprime'+mass+coup+'_'+sys
            #print 'hadd ' + outdir + 'Wprime'+mass+coup+'_'+sys+'.root ' + indirec +'/*' +'root'
            os.system('hadd '+outdir+'Wprime'+mass+coup+'_'+sys+'.root '+indirec+'/*.root')
