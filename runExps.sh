#! /bin/bash

MASTERURL = $(uname -n)

for FILE in ~/permanent-hdfs/datasets/*k/*
	do
		echo Running experiment for dataset ${FILE}
		./sparkR ExpWrapper.R $MASTERURL 1 10 $FILE 10 150 10 ${FILE}.results
	done
