#!/bin/bash
for IT in {11..50}
do
for DIM in 4 # 4000 8000
	do
		for REV in 90 # 85 80
		do
			Rscript genMat.R ${DIM}000 10 0.01 0.$REV  outputfile${DIM}000${REV}_${IT}
			# scp -i dfc.pem outputfile${DIM}000${REV}_${IT}_masked* root@ec2-54-201-141-204.us-west-2.compute.amazonaws.com:~/persistent-hdfs/datasets/${DIM}k/
		done
done
done
