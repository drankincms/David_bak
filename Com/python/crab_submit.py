#!/usr/bin/python

import os
import sys
import re
import fileinput
import datetime
import subprocess

import environment
from crabCommandLineParser import *

# banner
def banner():
    print '+-----------------------------------------------------'
    print '|'
    print '| Collect info about environment'
    print '| Submit grid CMSSW jobs with CRAB'
    print '|'
    print '| Uses: crab.cfg.template'
    print '|'
    print '| author: Gena Kukartsev'
    print '|'
    print '| (c) 2010'
    print '|'
    print '+-----------------------------------------------------'



def check_par_validity(cl_options):
    valid_ = True
    if (cl_options.dataset_name == None):
        valid_ = False

    return valid_

# globals
legend = '[CRAB submit]:'
package = 'LJMet/Com/'
crab_cfg_template = 'crab/crab.cfg.template'

############################################################
#
# main function
#

banner()

username = os.environ['USER']

# command line argument parser
options = parse_command_line(username)

if ( check_par_validity(options)==False ):
    print 'Invalid options, exiting...'
    #parser.print_help()
    sys.exit()

# collect and save the environment
_env = environment.environment(options)

    


# remove trailing data tier names from the dataset name
ds_truncated = options.dataset_name[0:options.dataset_name.rfind('/')]
print 'Dataset, trailing data tier names removed: ', ds_truncated

# Substitute slashes and dashes with underscores
base_name = ds_truncated.strip('/')
base_name = base_name.replace('/','_')
base_name = base_name.replace('-','_')
print 'Directory and file name base: ', base_name

# date and time
submit_date = datetime.datetime.today()
date_suffix = submit_date.strftime('%d%b%Y_%H%M%S')
print 'Date_time: ', date_suffix

# crab directory
crab_dir = options.crab_path+base_name+'/'+date_suffix
os.system('rm -rf '+crab_dir)
#os.system('mkdir -p '+crab_dir)
print 'Crab directory: ', crab_dir


print legend, 'CRAB dir: ', crab_dir

crab_cfg_template_file = open(_env.cmssw_base+'/src/'+package+crab_cfg_template, 'r')
crab_cfg_file = open(_env.cmssw_base+'/src/'+package+'test/crab_'+base_name+'_'+date_suffix+'_cfg.crab',"w")
for line in crab_cfg_template_file:
    #line=line.replace('DIRECTORY',dir)
    #line=line.replace('PREFIX',prefix)
    #line=line.replace('JOBID',str(i))
    crab_cfg_file.write(line)
crab_cfg_file.close()
crab_cfg_template_file.close()

#    if (options.info == False):
#        os.system('chmod u+x '+dir+'/'+prefix+'_'+str(i)+'.csh')
#        os.system('crab_submit '+dir+'/'+prefix+'_'+str(i)+'.crab')

    
# CRAB directory
print 'CRAB directory: ', crab_dir
 
    
    
_env.summary()
