#! /usr/bin/python 

import csv
import sys
import os.path


# Time #slices error iterations m n #revealed
# slices, iterations, RMSE, overhead, subproblem, time projection, time visible

INF = float("inf")
# Tuples: procs, iters, error, overhead, subprob time, proj time
SLICES = 0
ITERS = 1
ERR = 2
OH = 3
T1 = 4
T2 = 5

timeBudgets = [100, 125, 150, 175, 200, 250, 300, 400, 700]

minErrDict = {}

for dataFile in os.listdir("."):
	if not ".out" in dataFile or not "movie" in dataFile:
		continue
	i = int(dataFile[17]) # Get the number of the matrix
	if not i in minErrDict:
		minErrDict[i] = (INF, INF, INF, INF, INF, INF, INF, INF, INF)
	f = open(dataFile, 'r')
	r = csv.reader(f, delimiter='\t', quotechar='|')
	errTup = minErrDict[i]
	try:
		for row in r:
			map(str.strip, row)
			if row[0] == 'slices':
				continue
			time = float(row[OH]) + float(row[T1]) +float(row[T2])
			for b in [t for t in timeBudgets if t >= time]:
				if errTup[timeBudgets.index(b)] > float(row[ERR]):
					T = [0]*len(errTup)
					for j in range(len(errTup)):
						if j == timeBudgets.index(b):
							T[j] = float(row[ERR])
						else:
							T[j] = errTup[j]
					errTup = tuple(T)
	except:
		print "Probrem"
		continue
	minErrDict[i] = errTup
	f.close()

dataDict = {}
for i in range(1,8):	
	dataDictEntry = [0]*(1+len(timeBudgets))
	avTimeMult = 0
	avErrMult = 0
	for b in timeBudgets:
		f = open("movie_no"+str(i)+"_"+str(b)+".results", 'r')
		r = csv.reader(f, delimiter='\t', quotechar='|')
		for row in r:
			if row[0] == "slices":
				continue
			map(str.strip, row)
			time = float(row[OH]) + float(row[T1]) + float(row[T2])
			err = float(row[ERR])
			budgetInd = timeBudgets.index(b)
			# multiple of min err, multiple of budget
			dataDictEntry[budgetInd] = (err/minErrDict[i][budgetInd],time/b)
			avTimeMult += time/b
			avErrMult += err/minErrDict[i][budgetInd]
	dataDictEntry[-1] = (avErrMult/len(timeBudgets), avTimeMult/len(timeBudgets))
	print dataDictEntry[-1]
	dataDict[i] = tuple(dataDictEntry)

# Write out minima to file
o = open("movie_table.csv", 'w')
p = open("movie_averages.csv", 'w')
w = csv.writer(o, delimiter='\t', quotechar='|')
v = csv.writer(p, delimiter='\t', quotechar='|')
for i in range(1,8):
	w.writerow(list(minErrDict[i]))
	v.writerow(list(dataDict[i]))
o.close()
