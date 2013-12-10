#! /usr/bin/python

import csv
import sys
import os.path


# Time #slices error iterations m n #revealed
# slices, iterations, RMSE, overhead, subproblem, time projection, time visible

for i in range(1,6):
	o = open("4k_no"+str(i)+"_DataTable.csv", 'w')
	for dataFile in os.listdir("./results/gaussian/4k"):
		if "_"+str(i)+"_" in dataFile:
			continue
		f = open("./results/gaussian/4k/"+dataFile, 'rb')
		r = csv.reader(f, delimiter='\t', quotechar='|')
		try:
			for row in r:
				map(str.strip, row)
				if row[0] == 'slices':
					continue
				time = float(row[3]) + float(row[4]) +float(row[5])
#			for matFile in os.listdir("./datasets/movielens"):
#				if matFile in dataFile:
#					g = open("./datasets/movielens/"+matFile, 'rb')
#					g.readline()
#					S = g.readline().split(" ")
#					m = S[0]
#					n = S[1]
#					rev = float(S[2])/(int(m)*int(n))
#					g.close()
				m = dataFile[10:14]
				n = m
				rev = (100.0 - float(dataFile[14:16]))/100.0
				string = str(time)+"\t"+str(row[0])+"\t"+str(row[2])+"\t"+str(row[1])+"\t"+str(m)+"\t"+str(n)+"\t"+str(rev)+"\n"
				o.write(string)
		except:
			continue
		f.close()
	o.close()

	
