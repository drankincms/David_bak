#!/usr/bin/env python

import os
#import sys
#import re
#import fileinput
import subprocess
import fnmatch

from optparse import OptionParser
import datetime

from sets import Set
import glob

import ROOT



# banner
def banner():
    print '+-----------------------------------------------------'
    print '|'
    print '| data_sync.py'
    print '|'
    print '| Sync, merge and other actions with data files'
    print '|'
    print '| author: Gena Kukartsev'
    print '|'
    print '|'
    print '| Merge:'
    print '|        ../python/data_sync.py -m -i data/ -w \'*ljmet*root\'' 
    print '|'
    print '+-----------------------------------------------------'



banner()

#############################################
#
#   command line argument parser
#
add_help_option = "./data_sync.py -a ACTION -s SOURCE_DIR [other options]"

parser = OptionParser(add_help_option)

parser.add_option("-m", "--merge", dest="merge", default=False,
                  help="Merge root files", action="store_true")

parser.add_option("-s", "--sync", dest="sync", default=False,
                  help="Sync directories", action="store_true")

parser.add_option("-i", "--input-dir", dest="source_dir", default=None,
                  help="Input directory", action="store")

parser.add_option("-o", "--output-directory", dest="destination_dir", default=None,
                  help="Destination directory", action="store")

parser.add_option("-w", "--wildcard", dest="wildcard", default="*/all.root",
                  help="Wildcard for the files, with path", metavar="WILDCARD")

(options, args) = parser.parse_args()
#
##### end of command line parser


ROOT.gROOT.ProcessLine('.L ../root/hadd.C++')



def MergeFiles(file_list):
    #
    # zeroth item is 'hadd'
    # first item in the list is the output file name
    # (to stay compatible with hadd inputs)
    #

    _outfname = file_list[1]
    _infnames = file_list[2:len(file_list)]

    #print
    #print 'DEBUG'
    #print _outfname
    #print _infnames

    _outfile = ROOT.TFile(_outfname, 'recreate')

    _inflist = ROOT.TList()
    for _name in _infnames:
        _inflist.Add(ROOT.TFile.Open(_name))

    ROOT.MergeRootfile(_outfile, _inflist)

    return



def sync(options):
    # input validity check
    if (options.source_dir == None):
        print 'No valid source dir specified, cannot sync, exiting...'
        return
    if (options.destination_dir == None):
        print 'No valid destination dir specified, cannot sync, exiting...'
        return
    
    files = subprocess.Popen(['find', options.source_dir, '-wholename', options.wildcard], stdout=subprocess.PIPE).communicate()[0]
    file_list = files.splitlines()
    
    
    for file in file_list:
        tFile = file[len(options.source_dir):len(file)]
        tFile = tFile.strip('/')
        #print tFile

        dest_file = options.destination_dir.strip('/')+'/'+tFile

        if (os.path.isfile(dest_file)):
            print 'File', dest_file, 'exists'
        else:
            #print 'Need to copy ', file, 'to', dest_file
            dest_file_dir = os.path.dirname(dest_file)
            print 'Copying', file, 'to', dest_file
            #print 'Creating', dest_file_dir
            subprocess.Popen(['mkdir', '-p', dest_file_dir], stdout=subprocess.PIPE).communicate()[0]
            subprocess.Popen(['cp', file, dest_file], stdout=subprocess.PIPE).communicate()[0]


def merge(options):

    #print 'DEBUG', options
    
    # input validity check
    if (options.source_dir == None):
        print 'No valid source dir specified, nothing to merge, exiting...'
        return
    #files = subprocess.Popen(['find', options.source_dir, '-name', '*root'], stdout=subprocess.PIPE).communicate()[0]
    files = subprocess.Popen(['find', options.source_dir, '-wholename', options.wildcard], stdout=subprocess.PIPE).communicate()[0]
    file_list = files.splitlines()

    #print 'DEBUG', files
    #return

    # set of unique directories
    sDir = Set()
    
    for file in file_list:
        sDir.add( os.path.dirname(file) )

    #print sDir

    for dir in sDir:
        if (os.path.exists(dir+'/all.root')):
            print ''
            print 'Directory:', dir
            print 'Merged file all.root already exists, no merging left to do...'
            print ''
        else:
            print ''
            print 'Directory:', dir
            print 'Merging all root files in', dir, 'into all.root...'

            hadd_args=['hadd', dir+'/all.root']
            #_inFiles=glob.glob(dir+'/*.root')
            _inFiles=glob.glob(dir+'/'+options.wildcard)
            
            # sort input files, so the empty ones come at the end
            _size = []
            _index = 0
            for _file in _inFiles:
                _size.append((_index, os.path.getsize(_file)))
                _index += 1

            #print _size

            for _file in sorted(_size, key=lambda s: s[1], reverse=True):
                #print _file[0], _file[1]
                hadd_args.append(_inFiles[_file[0]])
            
            # print hadd_args
            # hadd is fucked up, better do merging ourselves
            #subprocess.Popen(hadd_args, stdout=subprocess.PIPE).communicate()[0]
            MergeFiles(hadd_args)


            print '...done'
            print ''
            

###########################
#
# Main function
#
print 'Source directory:', options.source_dir
print 'Destination directory:', options.destination_dir


if options.sync:
    sync(options)
elif options.merge:
    merge(options)
else:
    parser.print_help()
    
#
###### end of main function
