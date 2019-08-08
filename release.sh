#!/bin/bash

while true;
do
	condor_q -hold > hold.txt

	awk '{ print $1 }' hold.txt > hold2.txt

	awk '$1 ~ /\.0$/{print}' hold2.txt > hold3.txt

	while IFS= read -r line;
	do
		condor_qedit $line requestmemory $1;
		condor_qedit $line requestdisk $2;

	done < "hold3.txt"

	condor_release $USER
	sleep $3;

done;

