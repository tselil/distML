#!/bin/bash

for DIM in 1000 2000 4000 8000
do
	for REV in 90 85 80
	do
		Rscript genMat.R $DIM 10 0.01 0.$REV  outputfile${DIM}${REV}_2
		scp -i dfc.pem outputfile${DIM}${REV}_2* root@ec2-54-201-115-244.us-west-2.compute.amazonaws.com:~/persistent-hdfs/datasets/
	done
done
