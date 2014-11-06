#!//usr/bin/python

import subprocess

dbsql_query = '\'find dataset where dataset like /TTbar/Summer10-START36_V9_S09-v1/GEN-SIM-RECO\''

dbsql = subprocess.Popen(['dbsql', 'find', 'dataset'], stdout=subprocess.PIPE).communicate()[0]

print dbsql
