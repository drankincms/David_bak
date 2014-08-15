#!/usr/bin/python

import os
#import sys
#import re
import fileinput
import subprocess


# collect and save the environment setup
class environment:
    def __init__(self, options):
        self.options = options
        self.legend = '[Environment]:'
        self.username = os.environ['USER']
        self.pwd = os.environ['PWD']
        self.cmssw_base = os.environ['CMSSW_BASE']
        self.showtags = subprocess.Popen(['showtags'], stdout=subprocess.PIPE).communicate()[0]

        self.dbsql_query = '\'find dataset, dataset.release, dataset.tag, sum(block.numevents) where dataset like ' + self.options.dataset_name + '\''

        print self.legend, 'query DBS...'
        
        if options.debug > 0:
            print self.legend, self.dbsql_query

        self.dbsql = subprocess.Popen(['dbsql', 'find',
                                       'site,', 'dataset,',
                                       'dataset.release,', 'dataset.tag,',
                                       'sum(block.numevents)', 'where',
                                       'dataset', 'like',
                                       self.options.dataset_name], stdout=subprocess.PIPE).communicate()[0]
        self.sites = []
        self.dataset = None
        self.release = None
        self.tag = None
        self.nevents = None
        for line in self.dbsql.splitlines():
            splitline = line.split()

            if options.debug > 2:
                print self.legend, line

            if len(splitline)<2:
                continue
            
            if splitline[0]=='site' or splitline[0]=='Using':
                continue
            
            self.sites.append(splitline[0])
            self.dataset = splitline[1]
            self.release = splitline[2]
            self.tag = splitline[3]
            self.nevents = splitline[4]



    def summary(self, filename='stdout'):
        summary = []
        summary.append('')
        summary.append('===> Submission environment summary')
        summary.append('')
        summary.append('Submitted from '+self.pwd)
        summary.append('CMSSW local base dir: '+self.cmssw_base)
        #summary.append('CRAB directory: '+self.options.crab_dir_name)
        #summary.append('Output file prefix: '+options.dataset_prefix)
        summary.append('Dataset: '+self.options.dataset_name)

        summary.append('dbsql info:')
        summary.append('   Dataset: '+self.dataset)
        summary.append('   RECO release: '+self.release)
        summary.append('   RECO global tag: '+self.tag)
        summary.append('   Number of events: '+self.nevents)

        summary.append('Available on:')
        for site in self.sites:
            summary.append('   '+site)
        
        summary.append(self.showtags)
        summary.append('')

        summary_file_name = filename
        if (filename!='stdout' and filename!=None):
            summary_file = open(summary_file_name, 'w')
        for line in summary:
            if (filename=='stdout' or filename==None):
                print line
            else:
                summary_file.write(line+'\n')
        if (filename!='stdout' and filename!=None):
            summary_file.close()

        print self.legend, 'summary is written to', summary_file_name

    
 
    
    
