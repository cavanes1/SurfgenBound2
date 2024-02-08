# import modules
import numpy as np
from os import listdir
from os.path import isfile, join
import os

# read files
# existing names with old order
f = open("names.all", "r")
oldnamefile = f.readlines()
f.close()
# new order
f = open("names.new", "r")
newnamefile = f.readlines()
f.close()

# extract data
# old name list
oldnames = []
for line in oldnamefile:
    currdat = line.split()
    oldnames.append(currdat[-1])
# new name list
newnames = []
for line in newnamefile:
    newnames.append(line.split()[-1])

# identify relevant files
mypath = "./"
allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
grdfiles = []
for fl in allfiles:
    if "cartgrd." in fl:
        grdfiles.append(fl)
print("Relevant files:")
print(grdfiles)

# back up files
relfiles = grdfiles + ["energy.all", "geom.all"]
os.system("mkdir OLDDATA")
for fl in relfiles:
    os.system("mv " + fl + " OLDDATA")

# sort
for name in newnames:
    origpos = oldnames.index(name)
    for fl in grdfiles + ["geom.all"]:
        f = open("OLDDATA/" + fl, "r")
        startread = origpos*6
        writethis = f.readlines()[startread:startread + 6]
        f.close()
        f = open(fl, "a")
        for wt in writethis:
            f.write(wt)
        f.close()
    f = open("OLDDATA/energy.all", "r")
    writethis = f.readlines()[origpos]
    f.close()
    f = open("energy.all", "a")
    f.write(writethis)
    f.close()

os.system("mv names.all OLDDATA")
os.system("mv names.new names.all")

print("\nNOW YOU MUST RUN POINTS.PY\n")
