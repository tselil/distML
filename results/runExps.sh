#! /bin/bash

#ec2-54-201-152-242.us-west-2.compute.amazonaws.com


for I in 1 2
do
  for J in {2..10}
  do
    for FILE in ~/persistent-hdfs/datasets/${I}k/*90_${J}_masked.out
    do
      for SLICES in 2 4 6 8 10
      do
        echo Running $SLICES slice experiment for dataset ${FILE}
        ~/spark/sparkR ExpWrapper.R spark://ec2-54-201-150-247.us-west-2.compute.amazonaws.com:7077 $SLICES $SLICES $FILE 10 50 10 ${FILE}.slices${SLICES}.results
      done
    done
  done
done
