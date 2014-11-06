#!/usr/bin/env python
#
# LJMet assistant
#
# Module: CRAB manager action
#
# Author: Gena Kukartsev
#
# September 2010
#
############################################

import sys
import os
import datetime
import subprocess
import re

from assistAction import assistAction
from optparse import OptionGroup
from assistEnv import environment

class assistCrab(assistAction):

    def __init__(self, parser):
        self.legend = '[CRAB manager]:'
        self.name = 'crab'

        # location of the CRAB template and default python config
        self.package = 'LJMet/Com/'
        self.crab_cfg_template = 'crab/crab.cfg.template'
        #self.python_config_template = 'python/pat_cfg.py'
        self.python_config_template = 'python/pat_fromReco_384p2_cfg.py'

        # local disk space for CRAB dirs and such
        self.default_crab_path = '/uscms_data/d1/lpcljm/production/crab'

        # storage element defaults to FNAL
        self.default_storage_element = 'cmssrm.fnal.gov'
        
        # storage path base defaults to FNAL
        self.default_storage_path_base = '/srm/managerv2?SFN='
        
        self.default_user_remote_dir_base = 'PAT'

        self.default_events_per_job = 50000

        self.default_lumis_per_job = 25
        
        self.register()

        self.add_options(parser)
        


    #########################################################
    #
    # Configure the action
    #
    #
    def configure(self, options):
        print self.legend, 'configuring CRAB action...'

        # dataset specified?
        if options.dataset_name == None:

            print self.legend, 'dataset was not specified, exiting...'
            sys.exit(-1)

        else:
            print self.legend, 'using dataset', options.dataset_name


        # MC or data specified?
        if options.is_mc == None:

            print self.legend, 'Not specified whether data or Monte Carlo, exiting...'
            sys.exit(-1)

        else:
            if options.is_mc:
                print self.legend, 'dataset will be processed as Monte Carlo'
            else:
                print self.legend, 'dataset will be processed as data'

                # JSON file
                if options.json_file == None:
                    print self.legend, 'JSON file not specified, exiting...'
                    sys.exit(-1)
                else:
                    # check that the JSON file exists
                    if os.path.isfile(options.json_file):
                        json_file = options.json_file
                    else:
                        print self.legend, 'JSON file does not exist, exiting...'
                        sys.exit(-1)
                        



        # collect and save the environment
        self.env = environment(options)


        # CMSSW python config file specified?
        if options.python_config == None:

            python_config_template = self.env.cmssw_base+'/src/'+self.package+self.python_config_template
            print self.legend, 'using default python config template file - good'

        else:
            print self.legend, '*'
            print self.legend, '* WARNING!'
            print self.legend, '* you define python config template to be used as', options.python_config
            print self.legend, '* Are you sure that you do not want to use the default?'
            print self.legend, '* If you are not sure, drop -p option.'
            print self.legend, '* If you are sure, remember that a few options, like global tag,'
            print self.legend, '* are passed to the python config by CRAB via command line parameters'
            print self.legend, '* These option will be ignored unless you parse the command line'
            print self.legend, '* parameters in you python config'
            print self.legend, '*'
            
            python_config_template = options.python_config

        print self.legend, 'python config template file:', python_config_template


        #global tag
        if options.global_tag == None:
            print self.legend, 'no global specified, will use the tag used for the dataset RECO'
            global_tag = self.env.tag
        else:
            global_tag = options.global_tag

        print self.legend, 'global tag:', global_tag


        # remove trailing data tier names from the dataset name
        ds_truncated = options.dataset_name[0:options.dataset_name.rfind('/')]

        if options.debug > 0:
            print self.legend, 'dataset, trailing data tier names removed: ', ds_truncated


        # Substitute slashes and dashes with underscores
        self.base_name = ds_truncated.strip('/')
        self.base_name = self.base_name.replace('/','_')
        self.base_name = self.base_name.replace('-','_')

        if options.debug > 0:
            print self.legend, 'base name:', self.base_name


        # date and time
        submit_date = datetime.datetime.today()
        self.date_suffix = submit_date.strftime('%d%b%Y_%H%M%S')

        if options.debug > 0:
            print self.legend, 'date_time: ', self.date_suffix


        # crab directory
        if options.crab_path == None:
            print self.legend, 'no CRAB path specified, default will be used:'

            crab_path = self.default_crab_path
            
            print self.legend, 'CRAB path: ', crab_path

        else:
            print self.legend, 'Are you sure that it is a good idea to redefine'
            print self.legend, 'CRAB path? It is advisable to use default.'
            print self.legend, 'You have chosen CRAB path: ', crab_path
            crab_path = options.crab_path
            
        self.crab_dir_base = crab_path+'/'+self.base_name
        self.crab_dir = self.crab_dir_base+'/'+self.date_suffix

        print self.legend, 'CRAB directory: cleaning up and creating', self.crab_dir
        os.system('rm -rf '+self.crab_dir)
        os.system('mkdir -p '+self.crab_dir)
        os.system('chmod -R 775 '+self.crab_dir_base)


        # storage element
        if options.storage_element == None:

            storage_element = self.default_storage_element
            print self.legend, 'no storage element specified, will use default'
            print self.legend, 'USER.storage_element =', storage_element

        else:

            storage_element = options.storage_element
            

        # storage path
        if options.storage_path == None:

            storage_path = self.default_storage_path_base+'/resilient/'+self.env.username+'/'
            print self.legend, 'no storage path specified, will use default'
            print self.legend, 'USER.storage_path =', storage_path

        else:

            storage_path = self.default_storage_path_base+options.storage_path
            

        if options.total_njobs:
            total_njobs = options.total_njobs
            print self.legend, 'total number of jobs is specified to be', total_njobs
            print self.legend, 'it takes precedence over [events_per_job] and [lumis_per_job]'

        else:
            # events per job
            if options.events_per_job == None:

                events_per_job = self.default_events_per_job
                print self.legend, 'number of events per job is not specified, default will be used'
                print self.legend, 'CMSSW.events_per_job =', events_per_job
                print self.legend, '(only used for Monte Carlo, ignored for data)'
                
            else:
                
                events_per_job = options.events_per_job
            
            # lumis per job
            if options.lumis_per_job == None:
                
                lumis_per_job = self.default_lumis_per_job
                print self.legend, 'number of lumi sections per job is not specified, default will be used'
                print self.legend, 'CMSSW.lumis_per_job =', lumis_per_job
                print self.legend, '(only used for data, ignored for Monte Carlo)'
                
            else:
                
                lumis_per_job = options.lumis_per_job


        # total number of events
        if options.total_nevents == None:
            total_nevents = -1
        else:
            total_nevents = options.total_nevents
            

            

        # total number of lumis
        if options.total_nlumis == None:
            total_nlumis = -1
        else:
            total_nlumis = options.total_nlumis
            

        # user_remote_dir
        user_remote_dir = self.default_user_remote_dir_base+'/'+self.base_name+'/'+self.date_suffix

        
        
        # create the python config out of template
        # minimal additions are done to the file, only to ensure
        # safe transfer of the parameters
        with open(python_config_template) as python_config_template_file:
            
            # all files are created in the current dir
            python_config = self.env.pwd.rstrip('/')+'/'+self.base_name+'_'+self.date_suffix+'_cfg.py'
            print self.legend, 'creating the python config file', python_config
            python_config_file = open(python_config,"w")

            for line in python_config_template_file:
                if options.debug > 2:
                    print self.legend, line
                
                python_config_file.write(line)

            # override global tag
            python_config_file.write('''

####################################################
#
# automatically generated by LJMET assistant script
#

# define command line parameters
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ()
options.register ('globalTag',
                  None, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Global tag")
options.register ('data',
                  None, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,            # string, int, or float
                  "Data(1) or MC(0)")
# setup defaults
options.globalTag = None
# get and parse the command line arguments
options.parseArguments()
#
###### end of command line parameter parser ########

## global tag
process.GlobalTag.globaltag = cms.string(options.globalTag)

#
## end of LJMET assistant script additions
            ''')

            python_config_file.close()


        # create CRAB config
        template_filename = self.env.cmssw_base+'/src/'+self.package+self.crab_cfg_template

        if options.debug > 0:
            print self.legend, 'using CRAB config template file:', template_filename
            

        # open CRAB cfg template
        with open(template_filename) as crab_cfg_template_file:
            
            crab_cfg_file_name = 'crab_'+self.base_name+'_'+self.date_suffix+'.cfg'

            # all files are created in the current dir
            #crab_cfg_file = open(self.env.cmssw_base+'/src/'+package+'test/'+crab_cfg_file_name,"w")
            crab_cfg_file = open(self.env.pwd.rstrip('/')+'/'+crab_cfg_file_name,"w")


            # save file name - will need for job submission
            self.crab_cfg_file_name = crab_cfg_file_name


            
            for line in crab_cfg_template_file:

                keep_line = True
                
                is_at_fnal = False
                for site in self.env.sites:
                    if 'fnal.gov' in site:
                        is_at_fnal = True

                if options.use_server == None:
                    _useserver = 1
                else:
                    _useserver = options.use_server

                _scheduler = 'glite'

                if (is_at_fnal):
                    _scheduler = 'condor'
                    _useserver = 0
            
                line=line.replace('CONFIGNAME',crab_cfg_file_name)
                line=line.replace('USESERVER', str(_useserver))
                line=line.replace('SCHEDULER', _scheduler)
                line=line.replace('DATASET', self.env.dataset)
                line=line.replace('CMSSWCONFIG', python_config)
                line=line.replace('STORAGEELEMENT', storage_element)
                line=line.replace('STORAGEPATH', storage_path)
                line=line.replace('USERREMOTEDIR', user_remote_dir)
                line=line.replace('UIWORKINGDIR', self.crab_dir)
                line=line.replace('TOTALNUMBEROFEVENTS', str(total_nevents))
                line=line.replace('TOTALNUMBEROFLUMIS', str(total_nlumis))
                line=line.replace('GLOBALTAG', str(global_tag))

                # optional
                # number of jobs
                if options.total_njobs:
                    line=line.replace('NUMBEROFJOBS', str(options.total_njobs))
                    #line=line.replace('EVENTSPERJOB', '')
                    if 'events_per_job' in line:
                        keep_line = False
                    #line=line.replace('LUMISPERJOB', '')
                    if 'lumis_per_job' in line:
                        keep_line = False
                else:
                    #line=line.replace('NUMBEROFJOBS', '')
                    if 'number_of_jobs' in line:
                        keep_line = False
                    line=line.replace('EVENTSPERJOB', str(events_per_job))
                    line=line.replace('LUMISPERJOB', str(lumis_per_job))

                
                # has to be after number of jobs block above
                if options.is_mc:
                    #line=line.replace('#MONTECARLO#', '')
                    line=line.replace('ISDATA', '0')
                    #if '#DATA#' in line:
                    #    keep_line = False
                    if 'lumi_mask' in line:
                        keep_line = False
                    if 'total_number_of_lumis' in line:
                        keep_line = False
                    if 'lumis_per_job' in line:
                        keep_line = False
                else:
                    #line=line.replace('#DATA#', '')
                    line=line.replace('JSONFILE', json_file)
                    line=line.replace('ISDATA', '1')
                    if 'total_number_of_events' in line:
                        keep_line = False
                    if 'events_per_job' in line:
                        keep_line = False
                    
                if keep_line:
                    crab_cfg_file.write(line)

            crab_cfg_file.close()

        if options.debug > 2:
            self.env.summary()


            
            

        return {'status':'pass'}

    

    def execute(self, options):

        print self.legend, 'executing CRAB action...'

        os.system('crab -create -cfg '+self.crab_cfg_file_name)
        #create_log = subprocess.Popen(['crab', '-create',
        #                               '-cfg', self.crab_cfg_file_name],
        #                              stdout=subprocess.PIPE).communicate()[0]
        #if options.debug>0:
        #    print create_log


        # record bookkeeping info into the CRAB dir
        self.env.summary(self.crab_dir.rstrip('/')+'/summary.log')



        if options.testonly:
            print self.legend, 'test run only, will submit no jobs'
        else:
            os.system('crab -submit all -c '+self.crab_dir)
            
        os.system('ln -s '+self.crab_dir+' '+self.env.pwd.rstrip('/')+'/'+self.base_name+'_'+self.date_suffix)
        #print '''
        #os.system('crab -submit all -c '+self.crab_dir)
        #'''



        #jobs_created = re.search('Total of ([0-9]+) jobs created', create_log).group(1);


        #if jobs_created > 0:
        #    print self.legend, jobs_created, 'CRAB jobs created, submitting...'
        #    submit_log = subprocess.Popen(['crab', '-submit',
        #                                   '-c', self.crab_dir],
        #                                  stdout=subprocess.PIPE).communicate()[0]
        #    print submit_log
        #else:
        #    print self.legend, jobs_created, 'no CRAB jobs created, nothing to submit'



        return {'status':'pass'}


    def add_options(self, parser):
        
        group = OptionGroup(parser, "CRAB options",
                            '''
Usage: ./assist -a crab -d DATASETNAME [-p PYTHONCONFIG] [params]
                            ''')

# add num events and num lumis!!!!

        group.add_option("-d", "--dataset", dest="dataset_name", default=None,
                          help="Name of the dataset", metavar="DATASETNAME")
        
        group.add_option("--mc", dest="is_mc", default=None,
                         help="qualify dataset as Monte Carlo",
                         action="store_true", metavar="ISMC")
        
        group.add_option("--data", dest="is_mc", default=None,
                         help="qualify dataset as data",
                         action="store_false", metavar="ISMC")
        
        group.add_option("--test", dest="testonly", default=False,
                         help="Test run: create jobs but do not submit",
                         action="store_true", metavar="TESTONLY")
        
        group.add_option("-j", "--json-file", dest="json_file", default=None,
                          help="JSON file (list of good lumi sections)", metavar="JSONFILE")
        
        group.add_option("-p", "--python-config", dest="python_config", default=None,
                          help="CMSSW python config file", metavar="PYTHONCONFIG")
        
        group.add_option("-c", "--crab-path", dest="crab_path", default=None,
                          help="Base location of CRAB directories", metavar="CRABPATH")
        
        group.add_option("-r", "--storage-element", dest="storage_element", default=None,
                          help="CRAB storage element", metavar="STORAGEELEMENT")
        
        group.add_option("-s", "--storage-path", dest="storage_path", default=None,
                          help="CRAB storage path", metavar="STORAGEPATH")
        
        group.add_option("-t", "--global-tag", dest="global_tag", default=None,
                          help="Global tag", metavar="GLOBALTAG")
        
        group.add_option("-e", "--events-per-job", dest="events_per_job", default=None,
                          help="Number of events per job", metavar="EVENTSPERJOB")
        
        group.add_option("-n", "--total-number-of-events", dest="total_nevents", default=None,
                          help="Total number of events to process, -1 for all", metavar="TOTALNUMBEROFEVENTS")
        
        group.add_option("-N", "--total-number-of-lumis", dest="total_nlumis", default=None,
                          help="Total number of lumi sections to process, -1 for all", metavar="TOTALNUMBEROFLUMIS")
        
        group.add_option("-J", "--total-number-of-jobs", dest="total_njobs", default=None,
                          help="Total number of jobs to create", metavar="TOTALNUMBEROFJOBS")
        
        group.add_option("-l", "--lumis-per-job", dest="lumis_per_job", default=None,
                          help="Number of lumi sections per job", metavar="LUMISPERJOB")
        
        group.add_option("-S", "--use-server", dest="use_server", default=None,
                          help="Use crab server", metavar="USESERVER")

        #group.add_option("-t", action="store", help="Test option",
        #                 dest="test", default = None, metavar="TEST")

        parser.add_option_group(group)
        #parser.print_help()
