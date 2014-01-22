from pymol import cmd
#from glob import glob
#import build_seq
#import seq_convert

#files = glob.glob("/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/sequences/*")

#for file in glob("sekwencje/*"):
#file="AA"
#runbuildseq='run /home/wiktor/Zabawki/pymol_scripts/lib/build_seq.py'
#cmd.do(runbuildseq)
#cmd.do('run /home/wiktor/Zabawki/pymol_scripts/lib/build_seq.py')
#for aa in "DCAHWLGELVWCT": cmd._alt(string.lower(aa))
residues = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
for residue in residues:
     for i in range(3):
        cmd._alt(residue.lower())
#     cmd.set_geometry ??

#cmd.save(aa.pdb, "all")





#run /home/wiktor/Zabawki/pymol_scripts/lib/build_seq.py
#cmd.do("run seq_convert.py")
#builds='build_seq'+'AA'+',ss=helix'
#cmd.do(builds)
#pdb='AA.pdb'+','+'ala'
#cmd.save(pdb)
#cmd.delete("all")
#    cmd.load(file,'prot')
#    for a in cmd.index("elem s and bound_to elem s"):
#        if cmd.select("s1","%s`%d"%a) and \
#           cmd.select("s2","elem s and bound_to %s`%d"%a):
#            if cmd.select("(s1|s2) and not ?skip"):
#                cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
#                print '   ',round(cmd.dist("tmp","s1","s2"),3)
#                cmd.select("skip","s1|s2|?skip")

#cd Scripts
#run build_seq_phi_psi.py
#run build_seq.py
#from glob import glob
#cd phi-psi-sequences
#for a in glob("*"): cmd.do("build_seq"+' '+a), cmd.save("/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/structures/"+a+".pdb"), cmd.delete("all")

#from pymol import cmd,stored
#import build_seq
#import seq_convert
#import sys, os

#path=('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/tmprandom/')
#dirlst = listdir(path)

#for afile in dirlst:
#    fh = open(path+afile,'r')
#    for aline in fh.readlines():
#        if not aline startswith('>'):
#            build_seq(aline,None,None,None)
#	    cmd.save('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/'+ afile.pdb,(afile))



#cmd.do("run build_seq.py")
#fileopen=('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/tmprandom/randomA')

#for sequence in fileopen.readline():
#    for aa in sequence: cmd._alt(string.lower(aa))
    #cmd.save('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/Scripts/tmprandom/'+sequence+'.pdb',(selection))
   
#    cmd build_seq(aline)
#    cmd.save('/afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/'+afile.pdb, )

#for aa in "IYRSEH": cmd._alt(string.lower(aa))
#save /afs/pdc.kth.se/home/s/syazdi/Disc2/PROJECT/ile.pdb,(ile)
