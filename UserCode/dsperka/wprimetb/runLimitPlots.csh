#!/bin/sh

DIRS="./"
#echo $DIRS# x; y
-
title="e_mu+jets/BDT"

for sys in  allsys  
do
for coup in R #L RL
  do 
  for lepton in elec mu comb 
  #for lepton in comb 
    do

    title="e/mu+jets"
    if [[${lepton} == "mu"]]; then
     title="mu+jets/BDT"
    fi
    if [[${lepton}=="elec"]]; then
     title="e+jets/BDT"
    fi
    #FILES="${sys}_Wprime_${coup}_${lepton}"
    FILES="${lepton}"
    echo ${DIRS}  ${FILES} ${title}
    root -b -q plotLimits.C\(\"${DIRS}\",\"${FILES}\",\"${title}\"\)

done
done
done