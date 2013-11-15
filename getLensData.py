import sys

"""
if len(sys.argv < 2):
	print "python getLensData.py <input_filename> <output_filename>"
"""

outputfile = sys.argv[-1]
inputfile = sys.argv[-2]


maxRowEntry = 0
maxColEntry = 0
numNonzeros = 0
writeLines = []
f = open(inputfile, 'r')
for line in f:
	entryList = line.split('::')
	writeLines += [entryList]
	numNonzeros += 1
	maxRowEntry = max(int(entryList[0]), maxRowEntry)
	maxColEntry = max(int(entryList[1]), maxColEntry)
f.close()

movieIDs = [int(l[1]) for l in writeLines]
ids= (dict.fromkeys(movieIDs)).keys()
labels = range(len(ids))
pairs = zip(ids, labels)
idsToInts = dict(pairs)

o = open(outputfile, 'w')
o.write("%%MatrixMarket matrix coordinate real general\n")
o.write(str(maxRowEntry)+'\t'+str(len(idsToInts))+'\t'+str(numNonzeros)+'\n')
for l in writeLines:
	assert int(l[1]) in idsToInts
	o.write(l[0]+"\t"+str(idsToInts[int(l[1])])+'\t'+l[2]+'\n')
o.close()


