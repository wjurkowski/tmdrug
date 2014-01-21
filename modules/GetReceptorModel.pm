package GetReceptorModel;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(make_core_inputs make_core_single_inputs make_loops_inputs make_eval_inputs make_model_core make_model_core_m make_model_loops make_scoring select_core select_loops build_model);
$VERSION=1.0;

use lib '/home/wiktor/Zabawki/7TMpipe/TMpipe_mod';
use BaseFunctions;

sub build_model{
my ($ncore)=@_;
}

sub select_loops{
}

sub select_core{
my ($best,$minim);
my ($core_scorings,$scf)=@_;
for my $i (0..$#{$core_scorings}){
  for my $k (0..6){
	if($$scf==$k){
	  $best=0;
	  $minim=$$core_scorings[$k][0];
	  for my $i (1..4){
		if ($$core_scorings[$k][$i] < $minim) {
		$best=$i;
		$minim=$$core_scorings[$k][$i];}
	  }
	}
  }
}
return (\$best,\$minim);
}

sub make_scoring{
 my ($ncore)=@_;
 my (@core_scorings);
 #my (@molpdf, @DOPE, @ga341, @$i_bwDOPE, @i_bwmolpdf, @ZDOPE, @Zmolpdf);
 my $kb=8.62E-5;
 my $n=0;

 #retrieve scores
 my $make_core_out=$$ncore."_model_mult.log";
 open(MLOG, "<$make_core_out") or die "Can not open an input file: $!";
 seek(MLOG, 2, 6);
 1 while (defined($_ = <MLOG>));
 while (defined $_) {
    my @scory=split(/\s+/,$_);
	$molpdf[$n]=$scory[0];
	$DOPE[$n]=$scory[1];
	$ga341[$n]=$scory[2];
	$n++;
 }
 close MLOG;

 #construct weighted sum for score 
 my $avg_DOPE=avg(5,\@DOPE);
 my $avg_molpdf=avg(5,\@molpdf);
 my $stdev_DOPE=stdev(5,\@DOPE);
 my $stdev_molpdf=stdev(5,\@molpdf);
 my $sum_bwDOPE=0;
 my $sum_bwmolpdf=0;

 for my $i (0..4){
 $i_bwDOPE[$i]=exp(-$DOPE[$i]/($kb*273));
 $i_bwmolpdf[$i]=exp(-$molpdf[$i]/($kb*273));
 $sum_bwDOPE=$sum_bwDOPE+$i_bwDOPE[$i];
 $sum_bwmolpdf=$sum_bwmolpdf+$i_bwmolpdf[$i];
 $ZDOPE[$i]=($DOPE[$i]-$avg_DOPE)/$stdev_DOPE;
 $Zmolpdf[$i]=($molpdf[$i]-$avg_molpdf)/$stdev_molpdf;
 }

 for my $i (0..4){
 $core_scorings[0][$i]=$molpdf[$i];
 $core_scorings[1][$i]=$DOPE[$i];
 $core_scorings[2][$i]=$ga341[$i];
 $core_scorings[3][$i]=$Zmolpdf[$i]+$ZDOPE[$i]+$ga341[$i];
 $core_scorings[4][$i]=$i_bwDOPE[$i]/$sum_bwDOPE;
 $core_scorings[5][$i]=$i_bwmolpdf[$i]/$sum_bwmolpdf;
 $core_scorings[6][$i]=$ga341[$i]+$core_scorings[4][$i]+$core_scorings[5][$i];
 }
 return (@core_scorings);
}

sub make_model_loops{#consecutive output as inputs
my ($ncore)=@_;
my $loop_in6=$$ncore."_loop_refine.py";
my $loop_in4=$$ncore."_model_energies.py";

`mod9v8 $loop_in6`; 	#model loops; output: 
`mod9v8 $loop_in4`; 	#evaluate model; input: Target.B9999000*.pdb; output: Target.B9999000*.profile
}

sub make_model_core{#consecutive output as inputs
my ($ncore)=@_;
my $core_in1=$$ncore."_align.py";
my $core_in3=$$ncore."_model_mult.py";
my $core_in4=$$ncore."_model_energies.py";
#my $core_in5=$$ncore."_plot_profiles.py";

`mod9v8 $core_in1`; 	#align templates; output: fm00495.ali
`mod9v8 $core_in3`; 	#build model; input: Target-mult.ali; output: Target.B9999000*.pdb
`mod9v8 $core_in4`; 	#evaluate model; input: Target.B9999000*.pdb; output: Target.B9999000*.profile
#`mod9v8 $core_in5`; 	#print profiles; input: Target.B9999000*.profile; output: dope_profile.png
}

sub make_model_core_m{#consecutive output as inputs
my ($ncore)=@_;
my $core_in1=$$ncore."_salign.py";
my $core_in2=$$ncore."_align2d_mult.py";
my $core_in3=$$ncore."_model_mult.py";
my $core_in4=$$ncore."_model_energies.py";
#my $core_in5=$$ncore."_plot_profiles.py";

`mod9v8 $core_in1`; 	#align templates; output: fm00495.ali
`mod9v8 $core_in2`; 	#seq-structure alignemnt 
`mod9v8 $core_in3`; 	#build model; input: Target-mult.ali; output: Target.B9999000*.pdb
`mod9v8 $core_in4`; 	#evaluate model; input: Target.B9999000*.pdb; output: Target.B9999000*.profile
#`mod9v8 $core_in5`; 	#print profiles; input: Target.B9999000*.profile; output: dope_profile.png
}


sub make_loops_inputs{
my ($ncore,$normal_loop,$dope_loop,$ss,$receptor_model,$lewy,$prawy,$no_models,$cys1,$cys2)=@_;
my $loop_in6=$$ncore."_loop_refine.py";

#loop modelling
open(INP1, ">$loop_in6") if $loop_in6;
printf INP1 "\# Loop refinement of an existing model\n";
printf INP1 "from modeller import \*\n";
printf INP1 "from modeller.automodel import \*\n";
printf INP1 "\n";
printf INP1 "log.verbose\(\)\n";
printf INP1 "env = environ\(\)\n";
printf INP1 "\n";
printf INP1 "\# directories for input atom files\n";
printf INP1 "env.io.atom_files_directory = \'.\'\n";
printf INP1 "\n";
printf INP1 "\# Create a new class based on \'loopmodel\' so that we can redefine\n";
printf INP1 "\# select_loop_atoms \(necessary\)\n";
if ($$normal_loop == 1) {print INP1 "class myloop\(loopmodel\)\:\n";}
if ($$dope_loop == 1) {print INP1 "class myloop\(dope_loopmodel\)\:\n";}
printf INP1 "    \# This routine picks the residues to be refined by loop modeling\n";
printf INP1 "    def select_loop_atoms\(self\)\:\n";
printf INP1 "        \# 10 residue insertion \n";
printf INP1 "        return selection(self.residue_range\(\'$$lewy\', \'$$prawy\'\)\)\n";
printf INP1 "\n";
if ($$ss == 1) {
	printf INP1 "   def special_patches\(self, aln\)\: # A disulfide brigde between specified residues\n";
	printf INP1 "                 self.patch\(residue_type=\'DISU\', residues=\(self.residues\[\'$$cys1\'\], self.residues\[\'$$cys2\'\]\)\)\n";
}	
printf INP1 "m = myloop\(env, inimodel=\'$$receptor_model\', \# initial model of the target\n";
printf INP1 "%s%s%s\n","           sequence=\'loop.",$$lewy."-".$$prawy,"\',assess_methods=\(assess.DOPE, assess.GA341\)\)          # code of the target";
printf INP1 "\n";
printf INP1 "m.loop.starting_model= 1           \# index of the first loop model \n";
printf INP1 "m.loop.ending_model  = $$no_models          # index of the last loop model\n";
printf INP1 "m.loop.md_level = refine.very_fast  \# loop refinement method\n";
printf INP1 "\n";
printf INP1 "m.make\(\)\n";
printf INP1 "\n";
printf INP1 "\n";
close (INP1);
}

sub make_core_single_inputs{
my ($template,$chid,$ncore,$receptor_pir,$receptor_name)=@_;
my $core_in1=$$ncore."_align.py";
my $core_in3=$$ncore."_model.py";

#MINP1
open(INP1, ">$core_in1") if $core_in1;
print INP1 "from modeller import *\n";
print INP1 "\n"; 
print INP1 "log.verbose()\n";
print INP1 "env = environ()\n";
print INP1 "\n";
print INP1 "env.libs.topology.read(file=\'\$\(LIB\)\/top_heav.lib\'\)\n";
print INP1 "\n";
print INP1 "# Read aligned structure(s):\n";
print INP1 "aln = alignment(env)\n";
print INP1 "aln.append(file='fm00495.ali', align_codes='all')\n";
print INP1 "aln_block = len(aln)\n";
print INP1 "\n";
print INP1 "# Read aligned sequence(s):\n";
print INP1 "aln.append(file=\'$$receptor_pir\', align_codes=\'$$receptor_name\')\n";
print INP1 "\n";
print INP1 "# Structure sensitive variable gap penalty sequence-sequence alignment:\n";
print INP1 "aln.salign(output='', max_gap_length=20,\n";
print INP1 "           gap_function=True,   # to use structure-dependent gap penalty\n";
print INP1 "           alignment_type='PAIRWISE', align_block=aln_block,\n";
print INP1 "           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,\n";
print INP1 "           gap_penalties_1d=(-450, 0),\n";
print INP1 "           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),\n";
print INP1 "           similarity_flag=True)\n";
print INP1 "\n";
print INP1 "aln.write(file='Target-single.ali', alignment_format='PIR')\n";
print INP1 "#aln.write(file='Target-single.pap', alignment_format='PAP')\n";
close (INP1);

#MINP2
open(INP3, ">$core_in3") if $core_in3;
print INP3 "from modeller import *\n";
print INP3 "from modeller.automodel import *\n";
print INP3 "\n";
print INP3 "\n";
print INP3 "env = environ()\n";
print INP3 "a = automodel(env, alnfile='Target-single.ali',\n";
print INP3 "              knowns=(\'$$template\', \'$$chid]\'),sequence='Target',\n";
print INP3 "				assess_methods=(assess.DOPE, assess.GA341))\n";  
print INP3 "a.starting_model = 1\n";
print INP3 "a.ending_model = 5\n";
print INP3 "a.make()\n";
close (INP3);
}

sub make_core_inputs{
my ($templates,$templates_chid,$ncore,$receptor_pir,$receptor_name)=@_;
my @tempy=@{$templates};
my @chids=@{$templates_chid};
my $core_in1=$$ncore."_salign.py";
my $core_in2=$$ncore."_align2d_mult.py";
my $core_in3=$$ncore."_model_mult.py";

#MINP1
open(INP1, ">$core_in1") if $core_in1;
print INP1 "# Illustrates the SALIGN multiple structure/sequence alignment\n";
print INP1 "\n";
print INP1 "from modeller import *\n";
print INP1 "\n";
print INP1 "log.verbose()\n";
print INP1 "env = environ()\n";
print INP1 "env.io.atom_files_directory = './:../atom_files/'\n";
print INP1 "\n";
print INP1 "aln = alignment(env)\n";
print INP1 "for (code, chain) in (\n";
foreach my $i (0..$#tempy) {
	print MINP1 "(\'$tempy[$i]\', \'$chids[$i]\'),";
}
print INP1 "):\n";
print INP1 "    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))\n";
print INP1 "    aln.append_model(mdl, atom_files=code, align_codes=code+chain)\n";
print INP1 "\n";
print INP1 "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),\n";
print INP1 "                                    ((1., 0.5, 1., 1., 1., 0.), False, True),\n";
print INP1 "                                    ((1., 1., 1., 1., 1., 0.), True, False)):\n";
print INP1 "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,\n";
print INP1 "               rr_file='$(LIB)/as1.sim.mat', overhang=30,\n";
print INP1 "               gap_penalties_1d=(-450, -50),\n";
print INP1 "               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,\n";
print INP1 "               dendrogram_file='fm00495.tree',\n";
print INP1 "               alignment_type='tree', # If 'progresive', the tree is not\n";
print INP1 "                                      # computed and all structues will be\n";
print INP1 "                                      # aligned sequentially to the first\n";
print INP1 "               feature_weights=weights, # For a multiple sequence alignment only\n";
print INP1 "                                        # the first feature needs to be non-zero\n";
print INP1 "               improve_alignment=True, fit=True, write_fit=write_fit,\n";
print INP1 "               write_whole_pdb=whole, output='ALIGNMENT QUALITY')\n";
print INP1 "\n";
print INP1 "#aln.write(file='fm00495.pap', alignment_format='PAP')\n";
print INP1 "aln.write(file='fm00495.ali', alignment_format='PIR')\n";
print INP1 "\n";
print INP1 "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,\n";
print INP1 "           rr_file='$(LIB)/as1.sim.mat', overhang=30,\n";
print INP1 "           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),\n";
print INP1 "           gap_gap_score=0, gap_residue_score=0, dendrogram_file='fm00495.tree',\n";
print INP1 "           alignment_type='progressive', feature_weights=[0]*6,\n";
print INP1 "           improve_alignment=False, fit=False, write_fit=True,\n";
print INP1 "           write_whole_pdb=False, output='QUALITY')\n";
close (INP1);

#MINP2
open(INP2, ">$core_in2") if $core_in2;
print INP2 "from modeller import *\n";
print INP2 "\n"; 
print INP2 "log.verbose()\n";
print INP2 "env = environ()\n";
print INP2 "\n";
print INP2 "env.libs.topology.read(file='$(LIB)/top_heav.lib')\n";
print INP2 "\n";
print INP2 "# Read aligned structure(s):\n";
print INP2 "aln = alignment(env)\n";
print INP2 "aln.append(file='fm00495.ali', align_codes='all')\n";
print INP2 "aln_block = len(aln)\n";
print INP2 "\n";
print INP2 "# Read aligned sequence(s):\n";
print INP2 "aln.append(file=\'$$receptor_pir\', align_codes=\'$$receptor_name\')\n";
print INP2 "\n";
print INP2 "# Structure sensitive variable gap penalty sequence-sequence alignment:\n";
print INP2 "aln.salign(output='', max_gap_length=20,\n";
print INP2 "           gap_function=True,   # to use structure-dependent gap penalty\n";
print INP2 "           alignment_type='PAIRWISE', align_block=aln_block,\n";
print INP2 "           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,\n";
print INP2 "           gap_penalties_1d=(-450, 0),\n";
print INP2 "           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),\n";
print INP2 "           similarity_flag=True)\n";
print INP2 "\n";
print INP2 "aln.write(file='Target-multi.ali', alignment_format='PIR')\n";
print INP2 "#aln.write(file='Target-multi.pap', alignment_format='PAP')\n";
close (INP2);

#MINP3
open(INP3, ">$core_in3") if $core_in3;
print INP3 "from modeller import *\n";
print INP3 "from modeller.automodel import *\n";
print INP3 "\n";
print INP3 "\n";
print INP3 "env = environ()\n";
print INP3 "a = automodel(env, alnfile='Target-mult.ali',\n";
print INP3 "              knowns=\n";
foreach my $i ($#tempy) {
        print MINP1 "(\'$tempy[$i]\', \'$chids[$i]\'),";
}
print INP3 " sequence='Target',\n";
print INP3 "				assess_methods=(assess.DOPE, assess.GA341))\n";  
print INP3 "a.starting_model = 1\n";
print INP3 "a.ending_model = 5\n";
print INP3 "a.make()\n";
close (INP3);
}

sub make_eval_inputs{
my ($ncore,$no_models,$lewy,$prawy)=@_;
my $eval_in5=$$ncore."_plot_profiles.py";
my $dope_profile=$$ncore."_dope_profile.png";
my $eval_in4=$$ncore."_model_energies.py";

#energy calculation
open(INP4, ">$eval_in4") if $eval_in4;
#my $char='\%';
my $char  = chr(37);
printf INP4 "from modeller import *\n";
printf INP4 "from modeller.scripts import complete_pdb\n";
printf INP4 "\n";
printf INP4 "log.verbose\(\)    \# request verbose output\n";
printf INP4 "env = environ\(\)\n";
printf INP4 "env.libs.topology.read\(file=\'\$\(LIB\)\/top_heav.lib\'\) \# read topology\n";
printf INP4 "env.libs.parameters.read\(file=\'\$\(LIB\)\/par.lib\'\) \# read parameters\n";
printf INP4 "\n";
printf INP4 "%s%s%s\n","code = \"loop.",$$lewy."-".$$prawy,"\.IL00000001.pdb\"";
printf INP4 "mdl = complete_pdb\(env, code\)\n";
printf INP4 "s = selection\(mdl\)\n";
printf INP4 "s.assess_dope\(output=\'ENERGY_PROFILE NO_REPORT\', file=\'Target.profile\',\n";
printf INP4 "          normalize_profile=True, smoothing_window=15\)\n";
printf INP4 "for i in range\(1, $$no_models\)\:\n";
printf INP4 "    \# read model file\n";
printf INP4 "%s%s%s%s%s\n","    code = \"loop.",$$lewy."-".$$prawy,"\.BL",$char,"04d0001.pdb\" $char i";
printf INP4 "    mdl = complete_pdb\(env, code\)\n";
printf INP4 "    s = selection\(mdl\)\n";
printf INP4 "%s%s%s%s%s\n","    s.assess_dope\(output=\'ENERGY_PROFILE NO_REPORT\', file=\'loop.",$$lewy."-".$$prawy,"\.BL",$char,"04d0001.pdb.profile\' $char i ,";
printf INP4 "                  normalize_profile=True, smoothing_window=15\)\n";
close (INP4);

#MINP5
open(INP5, ">$eval_in5") if $eval_in5;
print INP5 "import pylab\n";
print INP5 "import modeller\n";
print INP5 "\n";
print INP5 "def get_profile(profile_file, seq):\n";
print INP5 "    \"\"\"Read `profile_file` into a Python array, and add gaps corresponding to\n";
print INP5 "       the alignment sequence `seq`.\"\"\"\n";
print INP5 "    # Read all non-comment and non-blank lines from the file:\n";
print INP5 "    f = file(profile_file)\n";
print INP5 "    vals = []\n";
print INP5 "    for line in f:\n";
print INP5 "        if not line.startswith('#') and len(line) > 10:\n";
print INP5 "            spl = line.split()\n";
print INP5 "            vals.append(float(spl[-1]))\n";
print INP5 "    # Insert gaps into the profile corresponding to those in seq:\n";
print INP5 "    for n, res in enumerate(seq.residues):\n";
print INP5 "        for gap in range(res.get_leading_gaps()):\n";
print INP5 "            vals.insert(n, None)\n";
print INP5 "    # Add a gap at position '0', so that we effectively count from 1:\n";
print INP5 "    vals.insert(0, None)\n";
print INP5 "    return vals\n";
print INP5 "\n";
print INP5 "e = modeller.environ()\n";
print INP5 "a = modeller.alignment(e, file='Target-mult.ali')\n";
print INP5 "\n";
print INP5 "#template = get_profile('1bdmA.profile', a['1bdmA'])\n";
print INP5 "model1 = get_profile('Target.B99990001.profile', a['Target'])\n";
print INP5 "model2 = get_profile('Target.B99990002.profile', a['Target'])\n";
print INP5 "model3 = get_profile('Target.B99990003.profile', a['Target'])\n";
print INP5 "model4 = get_profile('Target.B99990004.profile', a['Target'])\n";
print INP5 "model5 = get_profile('Target.B99990005.profile', a['Target'])\n";
print INP5 "\n";
print INP5 "# Plot the template and model profiles in the same plot for comparison:\n";
print INP5 "pylab.figure(1, figsize=(10,6))\n";
print INP5 "pylab.xlabel('Alignment position')\n";
print INP5 "pylab.ylabel('DOPE per-residue score')\n";
print INP5 "pylab.plot(model1, color='yellow', linewidth=2, label='Model1')\n";
print INP5 "pylab.plot(model2, color='orange', linewidth=2, label='Model2')\n";
print INP5 "pylab.plot(model3, color='red', linewidth=2, label='Model3')\n";
print INP5 "pylab.plot(model4, color='purple', linewidth=2, label='Model4')\n";
print INP5 "pylab.plot(model5, color='blue', linewidth=2, label='Model5')\n";
print INP5 "#pylab.plot(template, color='green', linewidth=2, label='Template')\n";
print INP5 "pylab.legend()\n";
print INP5 "pylab.savefig(\'$dope_profile\', dpi=65)\n";
close (INP5);
}
