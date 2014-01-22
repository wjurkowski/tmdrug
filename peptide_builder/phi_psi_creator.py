#This file takes a string/a sequence (output from sequence_permutator.py) and creates a new file contaning info about the phi and psi angles for each amino
#acid. The output files generated are saved in /afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/phi-psi-sequences/. The chosen phi and psi angles are
#-139 and 135, respectively.
import os.path
import sys
import re
from string import *

for a in range(1, len(sys.argv)):
	#file = open(sys.argv [a])
	file=open('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/sequences/' + str(a))
	writefile = open('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/phi-psi-sequences/'+str(a)+".txt", 'w')
	readfile = file.readlines()

	sequence=[]
	for a in readfile:
		sequence+=a.strip()
	for k in sequence:			
		if k == 'I':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)
		if k == 'Y':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)
		if k == 'R':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)
		if k == 'S':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)
		if k == 'E':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)
		if k == 'H':
			writefile.write(k+' '+'-139'+' '+'135'+"\n",)

