import csv
import sys
import numpy as np
from collections import defaultdict

### read case folder

caseFolder = sys.argv[1]
print ("Running testing in folder "+ caseFolder)
if caseFolder == "":
    print ("Please set case folder as a command line argument")

### input data

variables = ['T','p','U']
tEnd = 0.002
eps = 0.06

### needed functons ##########################

def readCSVFile (fileName):
    data = defaultdict(list)
    with open(fileName) as csvfile:
        rdr = csv.DictReader(csvfile)
        for row in rdr:
            for (k,v) in row.items():
                data[k].append(float(v))
    return data

def relativeError (ref, num):
    return np.fabs(num) if np.fabs(ref) < eps else np.fabs( (ref - num) / ref)

##############################################

print ("Error threshold: " + str(eps))

for var in variables:
    refDataFile = caseFolder + "/testing/ref/" + var + ".csv"
    numDataFile = caseFolder + "/postProcessing/sampleDict/" + str(tEnd) + "/"
    varNum = var
    if var == "U":
        numDataFile += "line_U.csv"
        varNum = "U_0"
    else:
        numDataFile += "line_T_p_rho.csv"

    refData = np.array(readCSVFile(refDataFile)[var])
    numData = np.array(readCSVFile(numDataFile)[varNum])

    ### compute minmax
    minRelError = relativeError(np.min(refData),np.min(numData))
    maxRelError = relativeError(np.max(refData),np.max(numData))
    print (var + "_min_error = " + str(minRelError))
    print (var + "_max_error = " + str(maxRelError))
    
    if (minRelError > eps or maxRelError > eps):
        print ("FAILED")
        sys.exit(1)

print ("OK")

    
### END
