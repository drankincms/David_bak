import os
import re
import glob
import sys

indir = '/uscms_data/d2/dsperka/8TeV/Samples/29May/'
outdir = '/uscms_data/d2/dsperka/8TeV/Samples/29May_All/'

list = [
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


if not os.path.exists(outdir):
    os.mkdir(outdir)

j = 0
for i in list:
   
    indirec = indir + i

    os.system('hadd ' + outdir + i +'.root ' + indirec +'/*' +'root')
 
    j += 1
    
