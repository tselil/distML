#! /bin/bash

SPARK_JAVA_OPTS+=" -Dspark.executor.memory=1g -Dspark.akka.frameSize=40"
export SPARK_JAVA_OPTS

MASTERURL=$(head -n 1 ~/spark-ec2/cluster-url)
echo $MASTERURL

for i in 1 2 3 4 6 7
do
    for TIME in 100 125 150 175 200 250 300 400 700
    do
        ./optimizer.py -m ~/persistent-hdfs/datasets/movielens/movielens10M_part${i}.mm -d dataTables/movie_no${i}_DataTable.csv -b 80000000 -t ${TIME} -u ${MASTERURL} -o optResults/movie_no${i}_${TIME}.results
    done
done
