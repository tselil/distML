import csv
import sys
import os.path


# Time #slices error iterations m n #revealed
# slices, iterations, RMSE, overhead, subproblem, time projection, time visible

o = open("dataTable.csv", 'w')
for dataFile in os.listdir("./results/"):
	f = open("./results/"+dataFile, 'rb')
	r = csv.reader(f, delimiter='\t', quotechar='|')
	try:
		for row in r:
			map(str.strip, row)
			if row[0] == 'slices':
				continue
			time = float(row[3]) + float(row[4]) +float(row[5])
			m = dataFile[10:14]
			n = m
			rev = 0.1
			o.write(str(time)+"\t"+str(row[0])+"\t"+str(row[2])+"\t"+str(m)+"\t"+str(n)+"\t"+str(rev)+"\n")
	except:
		continue
	f.close()
o.close()

	
