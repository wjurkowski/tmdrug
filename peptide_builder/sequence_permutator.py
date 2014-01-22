#This program takes a string and creates all possible permutations.
#The list consists of the amino acids I Y R S E H.
#The output files are saved in this following folder ('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/sequences/')
import sys, os, re
import random

def xselections(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss
#if __name__=="__main__":
#    for c in xselections(['I','Y','R','S','E','H'],3): print ''.join(c)

path=sys.argv[1]
i=1
if __name__=="__main__":
    for c in xselections(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],2):
	fileid='seq-'+repr(i)
	fh=open(path+fileid,'w')    	
	fh.write(''.join(c))	
	fh.close()
	i+=1	
	
