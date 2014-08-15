#!/usr/bin/env python

import os
import sys
import re
import fileinput

import datetime

import subprocess



legend = '[condor_submit]:'



# banner
def banner():
    print '''
+-----------------------------------------------------
|
| Submit FWLite LJMet/Com jobs to Condor
|
| authors: Gena Kukartsev, Michael Segala
|
| (c) 2010-2012
|
+-----------------------------------------------------
    '''
banner()



def RunDbsql(query, options):
    #
    # Equivalent of dbsql - runs a dbs query
    #

    # get the DBS instance URL
    _dbs_instance = 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
    if not options.dbs_instance:
        print legend, 'no DBS instance specified, default will be used'
    else:
        _dbs_instance = options.dbs_instance

    print legend, 'DBS instance:', _dbs_instance
        
    # command to run DBSQL queries
    dbs_cmd = os.environ['DBSCMD_HOME']+'/dbsCommandLine.py'


    print legend, 'quering DBS...'
    _pipe = subprocess.Popen(['python', dbs_cmd, '-c', 'search',
                              '--url', _dbs_instance, #'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet',
                              '--query', str(query)],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)

    stdout, stderr = _pipe.communicate()
    print legend, 'quering DBS... done'


    return stdout



# collect and save the environment setup
def GetEnvironment(options, dir):

    _environment = {}

    if len(options.dataset)<1:
        print legend, 'no dataset specified, exiting'
        sys.exit(-1)
        
    dataset = options.dataset[0]
        
    
    pwd = os.environ['PWD']

    cmssw_base = os.environ['CMSSW_BASE']

    showtags = subprocess.Popen(['showtags'], stdout=subprocess.PIPE).communicate()[0]

    summary = []
    summary.append('')
    summary.append('===> Submission environment summary')
    summary.append('')
    summary.append('Submitted from '+pwd)
    summary.append('CMSSW local base dir: '+cmssw_base)
    summary.append('Condor dir: '+dir)
    summary.append('Requested dataset: '+dataset)
    summary.append(showtags)
    summary.append('')

    summary_file_name = dir+'/summary.log'
    summary_file = open(summary_file_name, 'w')
    for line in summary:
        print line
        summary_file.write(line+'\n')
    summary_file.close()

    print 'Summary is written to', summary_file_name

    return _environment



def MakeCondorConfig(options, dir, prefix, jobid, input_files):
    #
    # create condor config
    #

    _config = '''
universe = vanilla
Executable = DIRECTORY/PREFIX_JOBID.csh
Requirements   =  OpSys == "LINUX" && (Arch =="INTEL" || Arch =="x86_64")
Should_Transfer_Files = YES
When_To_Transfer_Output = ON_EXIT
transfer_input_files = INPUT_FILES
Output = PREFIX_JOBID.stdout
Error = PREFIX_JOBID.stderr
Log = PREFIX_JOBID.condor.log

initialdir = INITIAL_DIR

Queue = 1
+UseSL5 = True
    '''

    _config = _config.replace('DIRECTORY', dir)
    _config = _config.replace('PREFIX', prefix)
    _config = _config.replace('JOBID', jobid)

    _inFiles = ''
    firstEntry = True
    for file in input_files:
        if not firstEntry:
            _inFiles += ', '
        _inFiles += file
        firstEntry = False

    _config = _config.replace('INPUT_FILES', _inFiles)
    _config = _config.replace('INITIAL_DIR', dir)
        

    return _config



def MakeCshellConfig(options, python_cfg):
    #
    # create csh config
    #

    if not options.application:
        print legend, 'no fwlite application specified, exiting'
        sys.exit(-1)
        

    _config = '''#!/bin/csh

set pwd = $cwd

source /uscmst1/prod/sw/cms/cshrc prod
setenv PATH ${{CMS_PATH}}/common:/bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog
echo ${{PATH}}

cd {cmssw_base}/src/

setenv SCRAM_ARCH {scram_arch}
eval `scramv1 runtime -csh`

#
#_____ check if we are running on a condor or on fbsng batch system _____
#
if (${{?_CONDOR_SCRATCH_DIR}}) then
  set BSCRATCH=${{_CONDOR_SCRATCH_DIR}}
  /bin/cat /etc/redhat-release
  echo "Batch system: Condor"
  echo "SCRAM_ARCH:   "${{SCRAM_ARCH}}
else if (${{?FBS_SCRATCH}}) then
  set BSCRATCH=$FBS_SCRATCH
  /bin/cat /etc/redhat-release
  echo "Batch system: FBSNG"
  echo "SCRAM_ARCH:   "${{SCRAM_ARCH}}
else
    echo "Unknown Batch System"
    exit
endif

cd $pwd

rehash
{fwlite} {python_config}
    '''

    _config = _config.format(cmssw_base=os.environ['CMSSW_BASE'],
                             scram_arch=os.environ['SCRAM_ARCH'],
                             fwlite=options.application,
                             python_config=os.path.basename(python_cfg))

    return _config



def MakePythonConfigCore(template, files, nevents, nskip, jobid):
    #
    # create python config from template
    #

    _config = ''

    for line in template.splitlines():

        # remove inline input definition
        if line.find('input_module')>-1:
            continue

        if line.find('process.inputs.nEvents')>-1:
            continue
            
        if line.find('process.inputs.skipEvents')>-1:
            continue
            
        if line.find('Input files')>-1:
            continue
            
        if line.find('btag_cond_file')>0:
            # remove btag conditions file since it will be transfered
            _line = line.split("=")
            if len(_line)>1:
                _config += '    btag_cond_file = cms.string(\''+os.path.basename(_line[1].strip().strip('\''))+'\n'
            continue
            
        if line.find('JsonFile')==0:
            # remove JSON file path since it will be transfered
            _line = line.split("=")
            if len(_line)>1:
                _config += 'JsonFile = \''+os.path.basename(_line[1].strip().strip('\''))+'\'\n'
            continue
            
        if line.find('outputName')>-1:
            _base = line.split('\'')
            if len(_base)>1:
                _file = _base[1] + '_' + jobid
                _config += '    outputName = cms.string(\''+_file+'\'),\n'
            continue

        _config += line + '\n'


    _input = '''
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes

process.inputs = cms.PSet (
    nEvents    = cms.int32({0}),
    skipEvents = cms.int32({1}),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring()
    )
    
'''

    for _file in files:
        _input += 'process.inputs.fileNames.extend(['

        # local or DBS file names
        if _file.find('pnfs')>0:
            # local file
            
            _input += '\'dcap://cmsgridftp.fnal.gov:24125'
            _newfname = _file.replace('pnfs','pnfs/fnal.gov/usr')
            _input += _newfname

        else:

            _input += '\'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11'
            _input += _file

        _input += '\'])\n\n'


    _input = _input.format( str(nevents), str(nskip) )


    _config += _input


    return _config



def MakePythonConfig(options, template, files, nevents, nskip, part_sum, jobid):
    #
    # create python config from template
    #

    _files = []

    import bisect

    # indeces of first and last needed files
    _iBegin = bisect.bisect_right(part_sum, nskip)
    _iEnd = min(bisect.bisect_left(part_sum, nskip+nevents), len(part_sum)-1)

    #print _iBegin, _iEnd 

    _added = 0
    _offset = -1  # events to skip from first file
    for _index in range(_iBegin, _iEnd+1):

        _files.append(files[_index][0])

        _added += files[_index][1]

    _before = 0 # events in skipped files
    if _iBegin > 0:
        _before = part_sum[_iBegin-1]

    # estimate the offset for first file in list
    _offset = nskip-_before


    return MakePythonConfigCore(template, _files, nevents, _offset, jobid)



def SubmitFromDbs(options):
    #
    # this function does everything
    #

    # query dataset name from options
    if len(options.dataset)<1:
        print legend, 'no dataset specified, exiting'
        sys.exit(-1)
        
    _dataset = options.dataset[0]


    # condor directory
    # specified_dir/dataset/date_version
    _dir = None
    if options.output_dir:
        _dir = options.output_dir
    else:
        _dir = './data'

    _dir = os.path.abspath(_dir)

    # check specified dir
    if os.path.exists(_dir):

        # create dataset subdir
        _dir += '/'+_dataset.strip('/').replace('/','_').replace('-','_')
        if not os.path.exists(_dir):
            os.makedirs(_dir)
            print legend, _dir, 'created'
        else:
            print legend, _dir, 'exists, will create output dir there'

        #create output dir as date_version
        from datetime import datetime
        _date = datetime.now()
        _sDate = _date.strftime('%d%b%Y')
        _ver = 0
        while(True):
            _oDir = '/'+_sDate+'v'+str(_ver)
            if os.path.exists(_dir+_oDir):
                _ver += 1
                continue
            _dir = _dir + _oDir
            os.makedirs(_dir)
            break

    else:
        print legend, _dir, 'does not exist, cannot create output directory, exiting'
        sys.exit(-1)
        

    print legend, 'condor dir:', _dir


    # save and print some environment summary
    _env = GetEnvironment(options, _dir)
    


    # query DBS
    _query = 'find dataset, file.name, file.numevents where dataset = ' + _dataset
    #print legend, 'DBSQL query:', _query
    _dbsout = RunDbsql(_query, options)
    _dbs = _dbsout.strip().split('\n')
    _ds = None
    _multiple_ds = False
    _files = set()
    _file_nev = [] #for each [file, nev, partial_sum(nev)]
    _part_nev = [] #partial event sum
    _nev = 0
    for line in _dbs:

        if len(line)==0 or line[0] != '/':
            continue

        _cols = line.split()

        if _ds and _ds != _cols[0]:
            _multiple_ds = True
        _ds = _cols[0]

        _nev += int(_cols[2])

        _files.add(_cols[1])

        _file_nev.append([_cols[1], int(_cols[2]), int(_nev)])
        _part_nev.append( int(_nev) )


    #_nev = sum(_file_nev.values())

    if _multiple_ds:
        print legend, 'multiple datasets found, be more specific, exiting'
        sys.exit(-1)

    if not _ds:
        print legend, 'no files for dataset', _dataset, 'found, exiting'
        sys.exit(-1)

    print legend, 'dataset', _ds
    print legend, len(_files), 'files found'
    print legend, _nev, 'events total'



    # get python config template from file
    _py_template = ''
    if options.python_config:
        with open(options.python_config, 'r') as _pyFile:
            _py_template = _pyFile.read()


    # limit number of events if provided
    if options.nevents_total:
        _nev = min(_nev, int(options.nevents_total))
    

    print legend, _nev, 'events to be processed'


    # create condor files
    _nPerJob = 100000
    if options.nevents_job:
        _nPerJob = int(options.nevents_job)
    _nSkipped = 0
    _nJob = 0
    _prefix = 'ljmet'
    json_file = None
    input_files=[]
    for line in _py_template.splitlines():
        if line.find('JsonFile')==0:
            json_file = line.split('=')[1].strip().strip('\'')
    if json_file:
        subprocess.call(['cp', json_file, _dir+'/'])
        #input_files = [_dir+'/'+os.path.basename(json_file)]
        #input_files = [os.path.basename(json_file)]
        input_files.append(os.path.basename(json_file))

    # btag conditions file
    btag_file = None
    for line in _py_template.splitlines():
        if line.find('btag_cond_file')>0:
            #btag_file = line.split('=')[1].strip().strip('\'')
            btag_file = line.split('(')[1].strip(',').strip(')').strip().strip('\'')
    if btag_file:
        subprocess.call(['cp', btag_file, _dir+'/'])
        #input_files = [_dir+'/'+os.path.basename(json_file)]
        input_files.append(os.path.basename(btag_file))
    while _nSkipped < _nev:

        # python config
        _python_config = _dir+'/'+_prefix+'_{0}'.format(str(_nJob).zfill(4))+'_cfg.py'
        with open(_python_config, 'w') as _file:
            _file.write(MakePythonConfig(options,
                                         _py_template,
                                         _file_nev,
                                         _nPerJob,
                                         _nSkipped,
                                         _part_nev,
                                         str(_nJob).zfill(4)))

        # condor config
        input_files.append(os.path.basename(_python_config))
        _condor_config = _dir+'/ljmet_{0}'.format(str(_nJob).zfill(4))+'.condor'
        with open(_condor_config, 'w') as _file:
            _file.write(MakeCondorConfig(options, _dir, _prefix, str(_nJob).zfill(4), input_files))
        input_files.pop()

        # csh script
        _csh = _dir+'/ljmet_{0}'.format(str(_nJob).zfill(4))+'.csh'
        with open(_csh, 'w') as _file:
            _file.write(MakeCshellConfig(options, _python_config))
        subprocess.call(['chmod', 'a+x', _csh])


        # submit
        _cmd = 'condor_submit '+_condor_config
        if options.info:
            print legend, _cmd
        else:
            print _cmd.split()
            subprocess.call(_cmd.split())

        _nSkipped += _nPerJob
        _nJob += 1
    
        if options.info:
            print legend, _nJob, 'jobs would have been submitted'
        else:
            print legend, _nJob, 'jobs submitted'
    

    return



def SubmitFromDir(options):
    #
    # this function submits unpublished dataset
    #


    # still need dataset name for naming output etc.
    if len(options.dataset)<0:
        print legend, 'no dataset name specified, exiting'
        sys.exit(-1)

    _dataset = options.dataset[0]


    # if input dir specified, use its contents as input
    _inDir = options.input_dir[0]
    print legend, 'take input files from', _inDir
        


    # condor directory
    # specified_dir/dataset/date_version
    _dir = None
    if options.output_dir:
        _dir = options.output_dir
    else:
        _dir = './data'

    _dir = os.path.abspath(_dir)

    # check specified dir
    if os.path.exists(_dir):

        # create dataset subdir
        _dir += '/'+_dataset.strip('/').replace('/','_').replace('-','_')
        if not os.path.exists(_dir):
            os.makedirs(_dir)
            print legend, _dir, 'created'
        else:
            print legend, _dir, 'exists, will create output dir there'

        #create output dir as date_version
        from datetime import datetime
        _date = datetime.now()
        _sDate = _date.strftime('%d%b%Y')
        _ver = 0
        while(True):
            _oDir = '/'+_sDate+'v'+str(_ver)
            if os.path.exists(_dir+_oDir):
                _ver += 1
                continue
            _dir = _dir + _oDir
            os.makedirs(_dir)
            break

    else:
        print legend, _dir, 'does not exist, cannot create output directory, exiting'
        sys.exit(-1)
        

    print legend, 'condor dir:', _dir


    # save and print some environment summary
    _env = GetEnvironment(options, _dir)
    


    # look into input dir
    _filenames = os.listdir(_inDir)
    _filenames.sort()
    _files = set()
    _part_nev = [] #partial event sum
    _nev = 0

    #print 'DEBUG:', _files

    for _file in _filenames:

        # skip if not a root file
        if _file.split('.')[-1]!='root':
            continue

        _pathfile = _inDir.rstrip('/')+'/'+_file

        _files.add(_pathfile)


    #print legend, 'DEBUG!!!', _files
    



    # get python config template from file
    _py_template = ''
    if options.python_config:
        with open(options.python_config, 'r') as _pyFile:
            _py_template = _pyFile.read()


    # total number of events to process - mandatory
    if options.nevents_total:
        _nev = int(options.nevents_total)
    else:
        print legend, 'total number of events missing, exiting.'
        sys.exit(-1)
    

    print legend, _nev, 'events to be processed'


    # create condor files
    _nPerJob = 100000
    if options.nevents_job:
        _nPerJob = int(options.nevents_job)
    _nSkipped = 0
    _nJob = 0
    _prefix = 'ljmet'
    json_file = None
    input_files=[]
    for line in _py_template.splitlines():
        if line.find('JsonFile')==0:
            json_file = line.split('=')[1].strip().strip('\'')
    if json_file:
        subprocess.call(['cp', json_file, _dir+'/'])
        #input_files = [_dir+'/'+os.path.basename(json_file)]
        #input_files = [os.path.basename(json_file)]
        input_files.append(os.path.basename(json_file))

    # btag conditions file
    btag_file = None
    for line in _py_template.splitlines():
        if line.find('btag_cond_file')>0:
            #btag_file = line.split('=')[1].strip().strip('\'')
            btag_file = line.split('(')[1].strip(',').strip(')').strip().strip('\'')
    if btag_file:
        subprocess.call(['cp', btag_file, _dir+'/'])
        #input_files = [_dir+'/'+os.path.basename(json_file)]
        input_files.append(os.path.basename(btag_file))


    while _nSkipped < _nev:


        # python config
        _python_config = _dir+'/'+_prefix+'_{0}'.format(str(_nJob).zfill(4))+'_cfg.py'
        with open(_python_config, 'w') as _file:
            _file.write(MakePythonConfigCore(_py_template,
                                             _files,
                                             _nSkipped+_nPerJob,
                                             _nSkipped,
                                             str(_nJob).zfill(4)))
            #print legend, 'DEBUG:', 'will create python config here'

        # condor config
        input_files.append(os.path.basename(_python_config))
        _condor_config = _dir+'/ljmet_{0}'.format(str(_nJob).zfill(4))+'.condor'
        with open(_condor_config, 'w') as _file:
            _file.write(MakeCondorConfig(options, _dir, _prefix, str(_nJob).zfill(4), input_files))
        input_files.pop()

        # csh script
        _csh = _dir+'/ljmet_{0}'.format(str(_nJob).zfill(4))+'.csh'
        with open(_csh, 'w') as _file:
            _file.write(MakeCshellConfig(options, _python_config))
        subprocess.call(['chmod', 'a+x', _csh])


        # submit
        _cmd = 'condor_submit '+_condor_config
        if options.info:
            print legend, _cmd
        else:
            print _cmd.split()
            subprocess.call(_cmd.split())

        print _nSkipped, _nPerJob
        _nSkipped += _nPerJob
        _nJob += 1
    
    if options.info:
        print legend, _nJob, 'jobs would have been submitted'
    else:
        print legend, _nJob, 'jobs submitted'
    

    return







############################################################
#
# main function
#

# command line argument parser
from optparse import OptionParser
add_help_option = "./submit_condor.py [other options]"

parser = OptionParser()

parser.add_option("--output-dir", dest="output_dir",
                  help="Output directory")

parser.add_option("-p", "--prefix", dest="dataset_prefix", default='noprefix',
                  help="Prefix for the file names", metavar="DATASETPREFIX")

parser.add_option("-i", "--info", action="store_true", dest="info", default=False,
                  help="Don't submit anything, only print info")

parser.add_option("-c", "--python-config", dest="python_config", default=None,
                  action="store",
                  help="Python config file name")

parser.add_option("--nevents-job", dest="nevents_job", default='1000000',
                  action="store",
                  help="Events per job")

parser.add_option("--nevents-total", dest="nevents_total", default=None,
                  action="store",
                  help="Total number of events to process")

parser.add_option("--dbs-instance", dest="dbs_instance", default=None,
                  action="store",
                  help="DBS instance URL")

parser.add_option("--application", dest="application", default='ljmet',
                  action="store",
                  help="FWLite application to run")

# datasets
parser.add_option("-d", "--dataset", dest="dataset", default=[],
                  action="append",
                  help="Dataset names", metavar="INFILE")

# input directory (unpblished datasets)
parser.add_option("--input-dir", dest="input_dir", default=[],
                  action="append",
                  help="Input directories")


print legend, 'parsing command line options...',
(options, args) = parser.parse_args()
print 'done'



# proceed
if len(options.input_dir)>0:
    # submit from input dir (no DBS)
    SubmitFromDir(options)

elif len(options.dataset)>0:
    SubmitFromDbs(options)

else:
    print legend, 'neither input dir nor dataset name specified, exiting'
