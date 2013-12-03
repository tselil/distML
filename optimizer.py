import sys
import csv
from collections import defaultdict

# Price per slice
price = 1

# Hardcoded param flags
TIME_FLAG = 0
BUDGET_FLAG = 7
ERROR_FLAG = 2

# Factor to determine nearness of matrix sizes
sizeFactor = 2

# Check if an m1-by-n1 matrix with p1 entries revealed
# is "near" in size to an m2-by-n2 matrix with p2 entries revealed
def isNear(m1,n1,p1,m2,n2,p2):
	ratio = (m1*n1*p1)/(m2*n2*p2)
	return ratio >= 1/sizeFactor and ratio <= sizeFactor

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
	
	dataTuples = [tuple(map(float,row)+[float(row[1])*float(row[0])*price]) \
	                  for row in dataReader\
	                  if (float(row[0]) <= time and \
	                        float(row[1])*float(row[0])*price <= budget and \
	                        float(row[2]) <= error and \
	                        isNear(float(row[4]),float(row[5]),float(row[6]),m,n,p))]	
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

# Choose best params given request and optimizer data
def chooseParams(dataTuples, paramToMinimize, explore):
	config = 0
	if explore:
		probDict = getExploreProbs(dataTuples,paramToMinimize)
	else:
		config = min(dataTuples,key=lambda t: t[paramToMinimize])
	return config
		
# Update optimizer data after run

# Main - Run DFC with chosen params
