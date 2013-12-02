#!/bin/bash
for IT in {3..10}
do
for DIM in 1 2 # 4000 8000
	do
		for REV in 90 # 85 80
		do
			Rscript genMat.R ${DIM}000 10 0.01 0.$REV  outputfile${DIM}000${REV}_${IT}
			scp -i dfc.pem outputfile${DIM}000${REV}_${IT}_masked* root@ec2-54-201-100-85.us-west-2.compute.amazonaws.com:~/persistent-hdfs/datasets/${DIM}k/
		done
done
done
