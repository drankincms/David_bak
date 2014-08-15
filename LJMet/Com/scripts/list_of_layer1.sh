#!/bin/bash
#
# by Gena Kukartsev
#
# This script outputs a list of Layer1 root files, for which
# the full CRAB bookkeeping information is available
#
# Input parameter - path to the crab directory. Assumes only one output
# root file in SE per job
# 
# 1. Get all stdout files from the crab directory
# 2. Check if corresponding stderr exists
# 3. Grep the stderr file for the output location (assume /pnfs/cms/WAX prefix)
# 4. Check if the corresponding output root file exists, and output its path
#

if [ -z "$1" ]
then
    echo '';
    echo 'This script outputs a list of Layer1 root files, for which';
    echo 'the full CRAB bookkeeping information is available';
    echo '';
    echo 'Usage:';
    echo '        list_of_layer1.sh path_to_crab_directory';
    echo '';
    exit;
fi;

echo 'import FWCore.ParameterSet.Config as cms'
echo ''
echo 'def PATInput() : '
echo '    readFiles = cms.untracked.vstring()'

nfiles=0
tempfile=`mktemp`
prefix='ljmet'

for file_stdout in `find $1/res/ -iname 'CMSSW*stdout' | xargs ls`
do
    num=`echo $file_stdout  | sed 's/.*CMSSW_\(.*\)\.stdout/\1/'`

# check if the stderr file exists
    file_stderr=`echo $file_stdout | sed 's/\.stdout/\.stderr/'`
    if [ -e $file_stderr ]
    then

	#tail -1000 $file_stderr  >> $tempfile

#'dcap:///pnfs/fnal.gov/usr/cms/WAX/resilient/mike1886/TTbar-madgraph_Summer09_7Tev/MinimumBias_Commissioning10_GOODCOLL-v8/myDataPatTuple_1_1.root'
	#file_layer1_path=/pnfs/cms/WAX`grep 'export SE_PATH' $1/job/CMSSW.sh | sed 's/.*\(\/resilient.*\)/\1/'`
	file_layer1_path=/pnfs/cms/WAX`grep 'export SE_PATH' $1/job/CMSSW.sh | sed 's/.*\(\/11.*\)/\1/'`
	#file_layer1_path=/store/data`grep 'export SE_PATH' $1/job/CMSSW.sh | sed 's/.*\(\/store.*\)/\1/'`
	file_layer1=""

	file_layer1_disk="x"
	for ((crabnumber=0; crabnumber<10; crabnumber++));
	do
	  for file in $file_layer1_path"*_"$num"_"$crabnumber"_*.root"
	  do
	    if [ -e $file ]
	    then
		file_layer1=`echo $file  | sed 's/pnfs/pnfs\/fnal.gov\/usr/'`
		file_layer1_disk=$file
	    fi;
	  done

	done
	
	#echo $file_layer1_disk

	if [ -e $file_layer1_disk ]
	then
	    tail -1000 $file_stdout  >> $tempfile

	    if [ "$nfiles" -eq 0 ]
	    then
		echo '    readFiles.extend( ['
	    fi;
	    if [ "$nfiles" -ne 0 ]
	    then
		echo ','
	    fi;
	    if [ "$nfiles" -ne 254 ]
	    then
		echo -n '	'"'"'dcap://'$file_layer1"'"
		nfiles=$((nfiles+1))
	    fi;
	    if [ "$nfiles" -eq 254 ]
	    then
		echo ','
		echo '	'"'"'dcap://'$file_layer1"'"' ] )'
		nfiles=0
	    fi;
	fi;
    fi;
done;

if [ "$nfiles" -ne 254 ]
then
    echo ' ] )'
fi;


echo ''
echo '    return cms.Source("PoolSource",'
echo '                      debugVerbosity = cms.untracked.uint32(200),'
echo '                      debugFlag = cms.untracked.bool(True),'
echo '                      fileNames = readFiles'
echo '                      )'
echo '#'
echo '# The original number of events: '`cat $tempfile | grep 'TrigReport Events total' | sed 's/.*total\ =\ \([0-9]*\)\ .*/\1/'| (tr '\n' +; echo 0) | bc`
echo '# Number of passed events: '`cat $tempfile | grep 'TrigReport Events total' | sed 's/.*total\ =\ \([0-9]*\)\ passed\ =\ \([0-9]*\)\ .*/\2/' | (tr '\n' +; echo 0) | bc`
rm $tempfile
