#!/usr/bin/python

from optparse import OptionParser

#############################################
#
#   command line argument parser
#
def parse_command_line(username):
    add_help_option = "crab -d dataset"

    parser = OptionParser(add_help_option)

    parser.add_option("-d", "--dataset", dest="dataset_name", default=None,
                      help="Name of the dataset", metavar="DATASETNAME")

    parser.add_option("-p", "--python-config", dest="python_config", default=None,
                      help="CMSSW python config file", metavar="PYTHONCONFIG")

    parser.add_option("-c", "--crab-path", dest="crab_path", default='/uscms_data/d1/lpcljm/'+username+'/crab/',
                      help="Location of the CRAB directory", metavar="CRABPATH")

    parser.add_option("-s", "--storage-path", dest="storage_path", default='/resilient/'+username+'/',
                      help="CRAB storage path", metavar="STORAGEPATH")

    parser.add_option("-u", "--user-remote-dir-base", dest="user_remote_dir_base", default='PAT/',
                      help="CRAB user remote dir base", metavar="USERREMOTEDIRBASE")

    parser.add_option("-t", "--global-tag", dest="global_tag", default='default',
                      help="Global tag", metavar="GLOBALTAG")

    parser.add_option("-e", "--events-per-job", dest="events_per_job", default=25000,
                      help="Number of events per job", metavar="EVENTSPERJOB")

    (options, args) = parser.parse_args()

    return options
#
##### end of command line parser
