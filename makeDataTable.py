#! /usr/bin/python

import csv
import sys
import os.path


# Time #slices error iterations m n #revealed
# slices, iterations, RMSE, overhead, subproblem, time projection, time visible

o = open("movieDataTable.csv", 'w')
for dataFile in os.listdir("./results/movielens"):
	if "part5" in dataFile or "part6" in dataFile or "part7" in dataFile:
		continue
	f = open("./results/movielens/"+dataFile, 'rb')
	r = csv.reader(f, delimiter='\t', quotechar='|')
	try:
		for row in r:
			map(str.strip, row)
			if row[0] == 'slices':
				continue
			time = float(row[3]) + float(row[4]) +float(row[5])
			for matFile in os.listdir("./datasets/movielens"):
				if matFile in dataFile:
					g = open("./datasets/movielens/"+matFile, 'rb')
					g.readline()
					S = g.readline().split(" ")
					m = S[0]
					n = S[1]
					rev = S[2]
					g.close()
			#m = dataFile[10:14]
			#n = m
			#rev = 0.1
			o.write(str(time)+"\t"+str(row[0])+"\t"+str(row[2])+"\t"+str(row[1])+"\t"+str(m)+"\t"+str(n)+"\t"+str(rev))
	except:
		continue
	f.close()
o.close()

	
