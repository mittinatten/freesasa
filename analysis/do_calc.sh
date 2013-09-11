#!/bin/bash

# user parameters
cullpdb_list=cullpdb_pc20_res1.6_R0.25_d130730_chains2117
pdbdir=~/pdb/
bindir=../

function LR {
    for d in 5.0 2.5 1.0 0.5 0.3 0.2 0.1 0.01
    do
	outputdir=LR$t/d$d/
	if [ ! -d "$outputdir" ]; then mkdir $outputdir; fi
	echo
	echo "###### $outputdir"
	for P in `cat $cullpdb_list | cut -b1-4` 
	do 
	    p=`echo $P|awk '{print tolower($0)}'`
	    echo -n .
	    $bindir/calc_sasa -L -d $d $pdbdir/$p.pdb -t $t > $outputdir/$p.sasa
	done
	echo
   done
}

function SR {
    for n in 20 50 100 200 500 2000 5000
    do
	outputdir=SR$t/$n/
	if [ ! -d "$outputdir" ]; then mkdir $outputdir; fi
	echo
	echo "###### $outputdir"
	for P in `cat $cullpdb_list | cut -b1-4` 
	do 
	    p=`echo $P|awk '{print tolower($0)}'`
	    #echo -n "$p "; 
	    echo -n .
	    $bindir/calc_sasa -S -n $n $pdbdir/$p.pdb -t $t> $outputdir/$p.sasa 
	done 
	echo
    done 
}

t=2
LR
SR
t=1
LR
SR
