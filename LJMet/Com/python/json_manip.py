#!/usr/bin/env python
#########################################################################
#
# json_manip.py
#
# json_manip module
#
# Usage: ./analyze.py --get-json
#                     --dataset <dataset1> [--dataset dataset2]...
#                    [--in-file file1.json] [--in-file file2.json]...
#
#
# Authors:
#          Gena Kukartsev, July 2012
#
#
#########################################################################



import os
import subprocess
import sys
import json
import re
import copy
#import glob



legend = '[JSON]:'



global_debug = False


class JsonSetRep:


    def __init__(self):
        self.mMap = {}



    def AddJson(self, json):
        #
        # add standard JSON as in a CMS file
        #
        for run in json.keys():
            _lumis = set() 
            for interval in json[run]:
                for _lumi in range(interval[0], interval[1]+1):
                    _lumis.add(_lumi)
            self.mMap[run] = _lumis

        return



    def MakeJson(self):
        #
        # generate JSON in a standard CMS format
        #
        _json = {}

        _runs = self.mMap.keys()
        _runs.sort()
        for run in _runs:
            _jLumis = []
            _lLumi = list(self.mMap[run])
            _lLumi.sort()
            _jLumi = []
            _lastlumi = -1
            for _lumi in _lLumi:
                if _lumi==_lastlumi+1:
                    if len(_jLumi)<2:
                        _jLumi.append(_lumi)
                    else:
                        _jLumi[1]=_lumi
                else:
                    if len(_jLumi)==2:
                        _jLumis.append(_jLumi)
                        _jLumi = []
                    elif len(_jLumi)==1:
                        _jLumi.append(_jLumi[0])
                        _jLumis.append(_jLumi)
                        _jLumi = []
                    _jLumi.append(_lumi)

                _lastlumi=_lumi

            if len(_jLumi)<2:
                #print run, _jLumi
                _jLumi.append(_jLumi[0])

            _jLumis.append(_jLumi)
                
            _json[run] = _jLumis

        return _json


    def GetNRuns(self):
        return len(self.mMap.keys())



    def GetNLumis(self):
        _count = 0
        for run in self.mMap.keys():
            _count += len(self.mMap[run])

        return _count



    def AddRunLumi(self, run, lumi):
        if str(run) not in self.mMap.keys():
            self.mMap[str(run)] = set()

        self.mMap[str(run)].add(lumi)

        return



    def GetMap(self):
        return self.mMap


    
    def Union(self, js):
        other_map = js.GetMap()
        for run in other_map.keys():
            for lumi in other_map[run]:
                self.AddRunLumi(run, lumi)
        return
        


    def CheckRunLumi(self, run, lumi):

        if run not in self.mMap.keys():
            return False

        if lumi not in self.mMap[run]:
            return False

        return True


        
    def Intersection(self, js):
        
        _jsonSet = JsonSetRep()

        _runs_other = js.GetMap().keys()
        _runs_self  = self.mMap.keys()
        _runs = list(set(_runs_other+_runs_self))
        

        #debug
        print '[debug]:  self n runs:', self.GetNRuns()
        print '[debug]: other n runs:', js.GetNRuns()
        print '[debug]:   all n runs:', len(_runs)

        for run in _runs:
            if run not in _runs_self or run not in _runs_other:
                continue
            else:
                sLumi = copy.deepcopy(self.mMap[run])
                sLumi.update(js.GetMap()[run])
                for lumi in sLumi:
                    if self.CheckRunLumi(run, lumi) and js.CheckRunLumi(run, lumi):
                        _jsonSet.AddRunLumi(run, lumi)

        return _jsonSet
        


def RunDbsql(query):
    #
    # Equivalent of dbsql - runs a dbs query
    #

    dbs_cmd = os.environ['DBSCMD_HOME']+'/dbsCommandLine.py'

    _pipe = subprocess.Popen(['python', dbs_cmd, '-c', 'search',
                              '--url', 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet',
                              '--query', str(query)],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)

    stdout, stderr = _pipe.communicate()

    #print stdout

    #return stdout, stderr
    return stdout


    
def ReadJsonFile(filename):
    #
    # take JSON file name and read it in
    #
    
    print legend, 'reading JSON file', os.path.basename(filename), '...'
    _itext = str()
    with open(filename) as _ifile:
        for line in _ifile:
            _itext += line

    exec('_json = '+_itext)

    print legend, 'done'

    return _json


    
def ReadDbs_SingleQuery(dataset):
    #
    # read runs and lumis for a list of datasets from DBS
    # (attempt with a signle query for all datasets - chokes)
    #
    
    # find all matching datasets
    print legend, 'searching for all matching datasets...'
    _query = 'find dataset where'
    _condition = ''
    _firstentry = True
    for ds in dataset:
        if _firstentry:
            _condition += ' dataset like '+str(ds)
            _firstentry = False
        else:
            _condition += ' or dataset like '+str(ds)

    #print 'DEBUG', _query+_condition
    _dbsout = RunDbsql(_query+_condition)

    # split stdout in lines, drop empty lines at the beginning and end
    _datasets = _dbsout.strip().split('\n')

    # remove the output header
    while len(_datasets[0])==0 or _datasets[0][0] != '/':
        _datasets.pop(0)

    print legend, 'done'

    if len(_datasets)<1:
        print legend, 'No matching dataset found, cannot do anything'
        return

    if len(_datasets)>1:
        print legend, 'WARNING! More than one dataset found'
        print legend, 'WARNING! Are you sure you want to merge all of these?'

    print legend, 'matching datasets:'
    for ds in _datasets:
        print legend, '   '+ds


    # find run and lumi for the new JSON
    print legend, 'reading run and lumi from DAS...'
    _dbsout = RunDbsql('find run,lumi where'+_condition)    

    _lines = _dbsout.split('\n')
    #print 'DEBUG', _lines

    _json_set = JsonSetRep()
    _run = -1
    for _line in _lines:
        words = _line.split()

        # skip header and lines that are not result of query
        if len(words)!=2:
            continue
        try:
            int(words[0])
        except ValueError:
            continue


        _run = int(words[0].strip())
        _lumi = int(words[1].strip())
        _json_set.AddRunLumi(_run, _lumi)

    print legend, 'done'

    # debug
    #print _json_set
    return _json_set


    
def ReadDbs(dataset):
    #
    # read runs and lumis for a list of datasets from DBS
    #
    
    # find all matching datasets
    print legend, 'searching for all matching datasets...'
    _query = 'find dataset where'
    _condition = ''
    _firstentry = True
    for ds in dataset:
        if _firstentry:
            _condition += ' dataset like '+str(ds)
            _firstentry = False
        else:
            _condition += ' or dataset like '+str(ds)

    print 'DEBUG:', global_debug
    if global_debug:
        with open('data/json/devel/dbsql_ds6.txt', 'r') as _file:
            _dbsout = _file.read()
    else:
        # query datasets from DAS
        _dbsout = RunDbsql(_query+_condition)

    # debug: dump query output to file
    #with open('dbsql1.txt', 'w') as f1:
    #    f1.write(_dbsout)


    # split stdout in lines, drop empty lines at the beginning and end
    _datasets = _dbsout.strip().split('\n')

    # remove the output header
    while len(_datasets[0])==0 or _datasets[0][0] != '/':
        _datasets.pop(0)

    print legend, 'done'

    if len(_datasets)<1:
        print legend, 'No matching dataset found, cannot do anything'
        return

    if len(_datasets)>1:
        print legend, 'WARNING! More than one dataset found'
        print legend, 'WARNING! Are you sure you want to merge all of these?'

    print legend, 'matching datasets:'
    for ds in _datasets:
        print legend, '   '+ds


    # find run and lumi for the new JSON
    print legend, 'reading run and lumi from DAS...'

    _dbsout = None
    _firstentry = True
    for ds in dataset:
        print legend, '   querying', ds.strip()
        if global_debug:
            if _firstentry:
                _firstentry = False
                with open('data/json/devel/dbsql_ds6_rl.txt', 'r') as _file:
                    _dbsout = _file.read()
        else:
            if not _dbsout:
                _dbsout = RunDbsql('find run,lumi where dataset like '+ds)
            else:
                _dbsout += RunDbsql('find run,lumi where dataset like '+ds)

    # debug
    with open('dbsql2.txt', 'w') as f1:
        f1.write(_dbsout)

    _lines = _dbsout.split('\n')
    #print 'DEBUG', _lines

    _json_set = JsonSetRep()
    _run = -1
    for _line in _lines:
        words = _line.split()

        # skip header and lines that are not result of query
        if len(words)!=2:
            continue
        try:
            int(words[0])
        except ValueError:
            continue


        _run = int(words[0].strip())
        _lumi = int(words[1].strip())
        _json_set.AddRunLumi(_run, _lumi)

    print legend, 'done'

    # debug
    #print _json_set
    return _json_set


    
def GetJsonFile(options):
    #
    # main module function to be called from outside
    #

    print legend, 'preparing JSON file'

    global_debug = options.debug
    print global_debug

    # full list of runs and lumis
    result = {}

    # json from input JSON files
    #runlumi = {}
    jsonSetInput = JsonSetRep()

    # json from DBS
    #runlumi_dbs = {}
    #jsonSetDbs = JsonSetRep()

    # read parent JSON files if any
    if (options.in_file):
        for fname in options.in_file:
            jsonSetInput.AddJson( ReadJsonFile(fname) )

        #print json.dumps(jsonSetInput.MakeJson())

        _nruns = jsonSetInput.GetNRuns()
        _nlumis = jsonSetInput.GetNLumis()
        
        print legend, 'read', _nruns, 'unique runs from JSON files'
        print legend, 'read', _nlumis, 'unique lumi sections from JSON files'


    # find all matching datasets if specified
    if options.dataset:

        _datasets = []

        # try to read list of datasets from a file
        # (if only one dataset option specified, it might be a file name)
        if len(options.dataset)==1:
            try:
                print legend, 'trying to open a file with dataset names...'
                with open(options.dataset[0]) as infile:
                    for line in infile:
                        line.strip()
                        if line[0]!='/':
                            continue
                        _datasets.append(line)
                        
            except IOError as e:
                print legend, 'no file with dataset names found'
                _datasets = options.dataset

        else:
            _datasets = options.dataset
                
        #jsonSetDbs = ReadDbs(options.dataset)
        jsonSetDbs = ReadDbs(_datasets)
        
        #print json.dumps(jsonSetDbs.MakeJson())

        _nruns = jsonSetDbs.GetNRuns()
        _nlumis = jsonSetDbs.GetNLumis()

        print legend, 'read', _nruns, 'unique runs from DBS'
        print legend, 'read', _nlumis, 'unique lumi sections from DBS'


    while True:

        if options.dataset and options.in_file:
            jsonSet = jsonSetDbs.Intersection(jsonSetInput)
            break

        if options.dataset:
            jsonSet = jsonSetDbs
            break

        if options.in_file:
            jsonSet = jsonSetInput
            break

        print legend, 'no dataset or JSON file specified'
        break

    _nruns = jsonSet.GetNRuns()
    _nlumis = jsonSet.GetNLumis()
    
    print legend, _nruns, 'unique runs in the final JSON'
    print legend, _nlumis, 'unique lumi sections in the final JSON'
            
    _json = jsonSet.MakeJson()

    if options.out_file:

        print legend, 'writing the resulting JSON to file', options.out_file
        
        with open(options.out_file, 'w') as ofile:
            sys.stdout = ofile
            print json.dumps(_json)
            sys.stdout=sys.__stdout__
    else:
        print legend, json.dumps(_json)
        print legend, 'specify output file name if you want JSON saved'

    return
