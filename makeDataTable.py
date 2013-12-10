#! /usr/bin/python

import csv
import sys
import os.path


# Time #slices error iterations m n #revealed
# slices, iterations, RMSE, overhead, subproblem, time projection, time visible

for i in range(1,8):
	o = open("movie_no"+str(i)+"_DataTable.csv", 'w')
	for dataFile in os.listdir("./results/movielens"):
		if "part"+str(i) in dataFile:
			continue
		f = open("./results/movielens/"+dataFile, 'rb')
		r = csv.reader(f, delimiter='\t', quotechar='|')
		try:
			for matFile in os.listdir("./datasets/movielens"):
				if matFile in dataFile:
					g = open("./datasets/movielens/"+matFile, 'rb')
					g.readline()
					S = g.readline().split(" ")
					m = S[0]
					n = S[1]
					rev = float(S[2])/(int(m)*int(n))
					g.close()
			for row in r:
				print row
				map(str.strip, row)
				if row[0] == 'slices':
					continue
				if float(row[0]) > 9:
					print row[0]
					continue
				time = float(row[3]) + float(row[4]) +float(row[5])
				print time
#				m = dataFile[10:14]
#				n = m
#				rev = (100.0 - float(dataFile[14:16]))/100.0
				string = str(time)+"\t"+str(row[0])+"\t"+str(row[2])+"\t"+str(row[1])+"\t"+str(m)+"\t"+str(n)+"\t"+str(rev)+"\n"
				o.write(string)
		except:
			continue
		f.close()
	o.close()

	
