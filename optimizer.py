#! /usr/bin/python

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

#./optimizer.py -m ~/persistent-hdfs/datasets/2k/outputfile200090_2_masked.out -d dataTable.csv -b 800 -t 200 -u spark://ec2-54-201-131-211.us-west-2.compute.amazonaws.com:7707

import sys
from subprocess import call
import getopt
import csv
from collections import defaultdict
import random

# Price per slice
price = 1

# Hardcoded param flags
TIME_FLAG = 0
BUDGET_FLAG = 7
ERROR_FLAG = 2

# Global variables representing maxs and stuff
MAX_TIME = 60*60*24 # One day
MAX_ERROR = 2
MAX_BUDGET = 1000000*MAX_TIME # ONE MILLION DOLLARS

# Factor to determine nearness of matrix sizes
sizeFactor = 2

# Check if an m1-by-n1 matrix with p1 entries revealed
# is "near" in size to an m2-by-n2 matrix with p2 entries revealed
def isNear(m1,n1,p1,m2,n2,p2):
	ratio = (m1*n1*p1)/(m2*n2*p2)
	return ratio >= 1.0/sizeFactor and ratio <= sizeFactor

# Add the non-parameter values in two data tuples together
# Assume the two tuples have matching input params
def sumNonParams(a,b):
	return (a[0]+b[0],a[1],a[2]+b[2],a[3],a[4],a[5],a[6],a[7])
	 
# Average all the times and errors for tuples with the same input params
# i.e. the same (#slices, #iterations, rows, cols, revealed entries)
def averageTuples(dataTuples):
	# Make a dictionary of lists of data tuples with the same input params
	paramDict = defaultdict(list)	
	for row in dataTuples:
		rowParams = (row[1], row[3], row[4], row[5], row[6], row[7])
		paramDict[rowParams].append(row)
			
	 # Average all the tuples with the same params
	averagedTuples = []
	for paramKey in paramDict:
		l = len(paramDict[paramKey])
		s = reduce(sumNonParams,paramDict[paramKey])
		averagedTuples.append((s[0]/l,s[1],s[2]/l,s[3],s[4],s[5],s[6],s[7]))
	 
	return averagedTuples

# Read the relevant optimizer data
# we assume row = [ time, #slices, error, iterations, m, n, #revealed]
def loadData(fileName, time, budget, error, m, n, p):
	dataFile = open(fileName,'r')
	dataReader = csv.reader(dataFile,delimiter='\t',quotechar='|')
	dataTuples = [] 
	for row in dataReader:
		timeOK = float(row[0]) <= time
		budgetOK = float(row[1])*float(row[0])*price <= budget
		errorOK = float(row[2]) <= error
		sizeOK = isNear(float(row[4]),float(row[5]),float(row[6]),m,n,p)
		if (timeOK and budgetOK and errorOK and sizeOK):
			print row
			print [float(row[4]),float(row[5]),float(row[6]),m,n,p]
			dataTuples += [tuple(map(float,row)+[float(row[1])*float(row[0])*price])]
	dataFile.close()
	return dataTuples

# Takes a list of tuples and a param flag and returns
# a dictionary of (tuple,probability) pairs for exploration mode
def getExploreProbs(dataTuples, paramFlag):
	probDict = {}
	denom = sum([1.0/t[paramFlag] for t in dataTuples]) 
	for t in dataTuples: 
		probDict[t] = (1.0/t[paramFlag])/denom
	return probDict

# Sample a tuple from the distribution given by a dictionary
def sampleTuple(probDict):
	cutOff = random.random()
	total = 0.0
	for t in probDict:
		total += probDict[t]
		if total >= cutOff:
			 return t

# Choose best params given request and optimizer data
def chooseParams(dataTuples, paramToMinimize, explore):
	config = 0
	# In explore mode we sample a tuple of params with probability
	# based on the param to minimize.
	if explore:
		probDict = getExploreProbs(dataTuples,paramToMinimize)
		config = sampleTuple(probDict)	
	else:
		# If in exploit mode we just return the best tuple
		config = min(dataTuples,key=lambda t: t[paramToMinimize])
	return config
			
# Update optimizer data after run
def updateData(fileName,tup):
	dataFile = open(fileName,'a')
	dataWriter = csv.writer(dataFile,delimiter='\t',quotechar='|')
	dataWriter.writerow(list(tup[0:7]))
	dataFile.close()
	 
# Main - Run DFC with chosen params
NUM_PARAMS = 5
def main(argv):
	datafile = "dataTable.csv"
	matrixfile = "default"
	outputfile = ""
	time = "unset"
	error = "unset"
	budget = "unset"
	optParam = ERROR_FLAG
	explore = False
	masterURL = "local"
	try:
		opts, args = getopt.getopt(argv,"hd:m:o:b:t:e:xu:")
		if ("-b" in argv and "-t" in argv and "-e" in argv):
			   raise Exception('Too many optimization params')
		if ("-b" not in argv and "-t" not in argv and "-e" not in argv):
			   raise Exception('Not enough optimization params')
	except: 
		print './optimizer.py -m <matrix> -d <data> -b <max_budget> -t <max_time>\
													-e <max_error> -u <masterURL> (-x) (-o <outputfile>)'
		print 'Exactly two of -b,-t,-e must appear.\n'
		sys.exit(2)
	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-d"):
			dataFile = arg
			print dataFile
		elif opt in ("-m"):
			matrixFile = arg
			print matrixFile
		elif opt == "-o":
			outputfile = arg
		elif opt == "-b":
			budget = float(arg)
		elif opt == "-e":
			error = float(arg)
		elif opt == "-t":
			time = float(arg)
		elif opt == "-x":
			explore = True
		elif opt == "-u":
			masterURL = arg
	if outputfile == "":
		outputfile = matrixFile+'_out'
	
	if time == 'unset':
		optParam = TIME_FLAG 
	elif error == "unset":
		optParam = ERROR_FLAG
	elif budget == "unset":
		optParam = BUDGET_FLAG

	f = open(matrixFile,'r')
	f.readline()
	matInfo = f.readline().split(" ")
	f.close()
	m = int(matInfo[0])
	n = int(matInfo[1])
	p = float(matInfo[2])/(m*n)
	tupList = loadData(dataFile, time, budget, error, m, n, p)
	tupList = averageTuples(tupList)
	config = chooseParams(tupList, optParam, explore)
	
	 # Call DFC
	slices = config[1]
	iterations = config[3]
	cmd =["~/spark/sparkR", "~/distML/DFC.R", masterURL, str(slices), matrixFile, str(iterations), outputfile]
	print config
	print '\n\n\n\n\n'
	call(cmd, shell=True)
	 
if __name__=="__main__":
	main(sys.argv[1:])
