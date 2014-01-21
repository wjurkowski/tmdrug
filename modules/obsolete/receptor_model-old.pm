package receptor_model.pm

sub make_model{#consecutive output as inputs

my $modeller_inp1=$ncore."_salign.py";
my $modeller_inp2=$ncore."_align2d_mult.py";
my $modeller_inp3=$ncore."_model_mult.py";
my $modeller_inp4=$ncore."_evaluate_model.py";
my $modeller_inp5=$ncore."_plot_profiles.py";

`python $modeller_inp1`; 	#align templates; output: fm00495.ali
`python $modeller_inp2`; 	#structural alignment; input: Target.ali,fm00495.ali; output: Target-multi.ali
`python $modeller_inp3`; 	#build model; input: Target-mult.ali; output: Target.B9999000*.pdb
`python $modeller_inp4`; 	#evaluate model; input: Target.B9999000*.pdb; output: Target.B9999000*.profile
`python $modeller_inp5`; 	#print profiles; input: Target.B9999000*.profile; output: dope_profile.png
}

sub make_loops{#consecutive output as inputs

my $modeller_inp6=$ncore."_loop_refine.py";
my $modeller_inp7=$ncore."_model_energies.py";
my $modeller_inp8=$ncore."_evaluate_model.py";

`python $modeller_inp6`; 	#model loops; output: 
`python $modeller_inp7`; 	# input: ; output: 
`python $modeller_inp8`; 	#input: ; output: 
}

sub test_model{
}

sub make_receptor_inputs{
my ($template_pdb1,$ch_template_pdb1,$template_pdb2,$ch_template_pdb2,$ncore)=@_;
my $modeller_inp1=$ncore."_salign.py";
my $modeller_inp2=$ncore."_align2d_mult.py";
my $modeller_inp3=$ncore."_model_mult.py";
my $modeller_inp4=$ncore."_evaluate_model.py";
my $modeller_inp5=$ncore."_plot_profiles.py";
my $dope_profile=$ncore."_dope_profile.png";

#MINP1
open(MINP1, $modeller_inp1) if $modeller_inp1;
print MINP1 "# Illustrates the SALIGN multiple structure/sequence alignment";
print MINP1 "";
print MINP1 "from modeller import *";
print MINP1 "";
print MINP1 "log.verbose()"
print MINP1 "env = environ()"
print MINP1 "env.io.atom_files_directory = './:../atom_files/'"
print MINP1 "";
print MINP1 "aln = alignment(env)";
print MINP1 "for (code, chain) in ((\'$template_pdb1\', \'$ch_template_pdb1\'), (\'$template_pdb2\', \'$ch_template_pdb2\')):";
print MINP1 "    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))";
print MINP1 "    aln.append_model(mdl, atom_files=code, align_codes=code+chain)";
print MINP1 "";
print MINP1 "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),";
print MINP1 "                                    ((1., 0.5, 1., 1., 1., 0.), False, True),";
print MINP1 "                                    ((1., 1., 1., 1., 1., 0.), True, False)):";
print MINP1 "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,";
print MINP1 "               rr_file='$(LIB)/as1.sim.mat', overhang=30,";
print MINP1 "               gap_penalties_1d=(-450, -50),";
print MINP1 "               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,";
print MINP1 "               dendrogram_file='fm00495.tree',";
print MINP1 "               alignment_type='tree', # If 'progresive', the tree is not";
print MINP1 "                                      # computed and all structues will be";
print MINP1 "                                      # aligned sequentially to the first";
print MINP1 "               feature_weights=weights, # For a multiple sequence alignment only";
print MINP1 "                                        # the first feature needs to be non-zero";
print MINP1 "               improve_alignment=True, fit=True, write_fit=write_fit,";
print MINP1 "               write_whole_pdb=whole, output='ALIGNMENT QUALITY')";
print MINP1 "";
print MINP1 "#aln.write(file='fm00495.pap', alignment_format='PAP')";
print MINP1 "aln.write(file='fm00495.ali', alignment_format='PIR')";
print MINP1 "";
print MINP1 "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,";
print MINP1 "           rr_file='$(LIB)/as1.sim.mat', overhang=30,";
print MINP1 "           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),";
print MINP1 "           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',";
print MINP1 "           alignment_type='progressive', feature_weights=[0]*6,";
print MINP1 "           improve_alignment=False, fit=False, write_fit=True,";
print MINP1 "           write_whole_pdb=False, output='QUALITY')";

#MINP2
open(MINP2, $modeller_inp2) if $modeller_inp2;
print MINP2 "from modeller import *";
print MINP2 ""; 
print MINP2 "log.verbose()";
print MINP2 "env = environ()";
print MINP2 "";
print MINP2 "env.libs.topology.read(file='$(LIB)/top_heav.lib')";
print MINP2 "";
print MINP2 "# Read aligned structure(s):";
print MINP2 "aln = alignment(env)";
print MINP2 "aln.append(file='fm00495.ali', align_codes='all')";
print MINP2 "aln_block = len(aln)";
print MINP2 "";
print MINP2 "# Read aligned sequence(s):";
print MINP2 "aln.append(file='Target.ali', align_codes='Target')";
print MINP2 ""
print MINP2 "# Structure sensitive variable gap penalty sequence-sequence alignment:"
print MINP2 "aln.salign(output='', max_gap_length=20,";
print MINP2 "           gap_function=True,   # to use structure-dependent gap penalty";
print MINP2 "           alignment_type='PAIRWISE', align_block=aln_block,";
print MINP2 "           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,";
print MINP2 "           gap_penalties_1d=(-450, 0),";
print MINP2 "           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),";
print MINP2 "           similarity_flag=True)";
print MINP2 "";
print MINP2 "aln.write(file='Target-multi.ali', alignment_format='PIR')";
print MINP2 "#aln.write(file='Target-multi.pap', alignment_format='PAP')";

#MINP3
open(MINP3, $modeller_inp3) if $modeller_inp3;
print MINP3 "from modeller import *";
print MINP3 "from modeller.automodel import *"
print MINP3 ""
print MINP3 ""
print MINP3 "env = environ()";
print MINP3 "a = automodel(env, alnfile='Target-mult.ali',";
print MINP3 "              knowns=(\'$template_pdb1.$ch_template_pdb1\',\'$template_pdb2.$ch_template_pdb2\'), sequence='Target')";
print MINP3 "a.starting_model = 1";
print MINP3 "a.ending_model = 5";
print MINP3 "a.make()";

#MINP4
open(MINP4, $modeller_inp4) if $modeller_inp4;
print MINP4 "from modeller import *";
print MINP4 "from modeller.scripts import complete_pdb";
print MINP4 "";
print MINP4 "log.verbose()    # request verbose output";
print MINP4 "env = environ()";
print MINP4 "env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology";
print MINP4 "env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters";
print MINP4 "";
print MINP4 "# read model file";
print MINP4 "mdl1= complete_pdb(env, 'Target.B99990001.pdb')";
print MINP4 "mdl2= complete_pdb(env, 'Target.B99990002.pdb')";
print MINP4 "mdl3= complete_pdb(env, 'Target.B99990003.pdb')";
print MINP4 "mdl4= complete_pdb(env, 'Target.B99990004.pdb')";
print MINP4 "mdl5= complete_pdb(env, 'Target.B99990005.pdb')";
print MINP4 "";
print MINP4 "# Assess all atoms with DOPE:";
print MINP4 "s = selection(mdl1)"
print MINP4 "s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.B99990001.profile',"
print MINP4 "              normalize_profile=True, smoothing_window=15)"
print MINP4 "s = selection(mdl2)"
print MINP4 "s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.B99990002.profile',"
print MINP4 "              normalize_profile=True, smoothing_window=15)"
print MINP4 "s = selection(mdl3)"
print MINP4 "s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.B99990003.profile',"
print MINP4 "              normalize_profile=True, smoothing_window=15)"
print MINP4 "s = selection(mdl4)"
print MINP4 "s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.B99990004.profile',"
print MINP4 "              normalize_profile=True, smoothing_window=15)"
print MINP4 "s = selection(mdl5)"
print MINP4 "s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.B99990005.profile',"
print MINP4 "              normalize_profile=True, smoothing_window=15)"

#MINP5
open(MINP5, $modeller_inp5) if $modeller_inp5;
print MINP5 "import pylab";
print MINP5 "import modeller";
print MINP5 "";
print MINP5 "def get_profile(profile_file, seq):";
print MINP5 "    \"\"\"Read `profile_file` into a Python array, and add gaps corresponding to";
print MINP5 "       the alignment sequence `seq`.\"\"\"";
print MINP5 "    # Read all non-comment and non-blank lines from the file:";
print MINP5 "    f = file(profile_file)";
print MINP5 "    vals = []";
print MINP5 "    for line in f:";
print MINP5 "        if not line.startswith('#') and len(line) > 10:";
print MINP5 "            spl = line.split()";
print MINP5 "            vals.append(float(spl[-1]))";
print MINP5 "    # Insert gaps into the profile corresponding to those in seq:";
print MINP5 "    for n, res in enumerate(seq.residues):";
print MINP5 "        for gap in range(res.get_leading_gaps()):";
print MINP5 "            vals.insert(n, None)";
print MINP5 "    # Add a gap at position '0', so that we effectively count from 1:";
print MINP5 "    vals.insert(0, None)";
print MINP5 "    return vals";
print MINP5 "";
print MINP5 "e = modeller.environ()";
print MINP5 "a = modeller.alignment(e, file='Target-mult.ali')";
print MINP5 "";
print MINP5 "#template = get_profile('1bdmA.profile', a['1bdmA'])";
print MINP5 "model1 = get_profile('Target.B99990001.profile', a['Target'])";
print MINP5 "model2 = get_profile('Target.B99990002.profile', a['Target'])";
print MINP5 "model3 = get_profile('Target.B99990003.profile', a['Target'])";
print MINP5 "model4 = get_profile('Target.B99990004.profile', a['Target'])";
print MINP5 "model5 = get_profile('Target.B99990005.profile', a['Target'])";
print MINP5 "";
print MINP5 "# Plot the template and model profiles in the same plot for comparison:";
print MINP5 "pylab.figure(1, figsize=(10,6))";
print MINP5 "pylab.xlabel('Alignment position')";
print MINP5 "pylab.ylabel('DOPE per-residue score')";
print MINP5 "pylab.plot(model1, color='yellow', linewidth=2, label='Model1')";
print MINP5 "pylab.plot(model2, color='orange', linewidth=2, label='Model2')";
print MINP5 "pylab.plot(model3, color='red', linewidth=2, label='Model3')";
print MINP5 "pylab.plot(model4, color='purple', linewidth=2, label='Model4')";
print MINP5 "pylab.plot(model5, color='blue', linewidth=2, label='Model5')";
print MINP5 "#pylab.plot(template, color='green', linewidth=2, label='Template')";
print MINP5 "pylab.legend()";
print MINP5 "pylab.savefig(\'$dope_profile\', dpi=65)";
}

sub make_loops_inputs{
my $modeller_inp6=$ncore."_loop_refine.py";
my $modeller_inp7=$ncore."_model_energies.py";
my $modeller_inp8=$ncore."_evaluate_model.py";
my ($ncore)=@_;

#loop modelling
open(LINP1, $modeller_inp6) if $modeller_inp6;
print LINP1 �\# Loop refinement of an existing model�;
print LINP1 �from modeller import *�;
print LINP1 �from modeller.automodel import *�;
print LINP1 ��;
print LINP1 �log.verbose()�;
print LINP1 �env = environ()�;
print LINP1 ��;
print LINP1 �\# directories for input atom files�;
print LINP1 �env.io.atom_files_directory = './:../atom_files'�;
print LINP1 ��;
print LINP1 �\# Create a new class based on 'loopmodel' so that we can redefine�;
print LINP1 �\# select_loop_atoms (necessary)�;
if ($normal_loop=1) print LINP1 �class myloop(loopmodel):�;
if ($dope_loop=1) print LINP1 �class myloop(dope_loopmodel):�;
print LINP1 �    # This routine picks the residues to be refined by loop modeling�;
print LINP1 �    def select_loop_atoms(self):�;
print LINP1 �        # 10 residue insertion �;
print LINP1 �        return selection(self.residue_range('141', '175'))�;
print LINP1 ��;
if ($disulfids=1) {
	print LINP1 �   def special_patches(self, aln): # A disulfide brigde between residues 15 and 187�;
	print LINP1 �                 self.patch(residue_type='DISU', residues=(self.residues['71'], self.residues['159']))�;
}	
print LINP1 �m = myloop(env, inimodel='multi2-Target.B99990001.pdb', \# initial model of the target�;
print LINP1 �           sequence='Target')          # code of the target�;
print LINP1 ��;
print LINP1 �m.loop.starting_model= 1           # index of the first loop model �;
print LINP1 �m.loop.ending_model  = 200          # index of the last loop model�;
print LINP1 �m.loop.md_level = refine.very_fast  \# loop refinement method�;
print LINP1 ��;
print LINP1 �m.make()�;
print LINP1 ��;
print LINP1 ��;

#energy calculation
open(LINP2, $modeller_inp7) if $modeller_inp7;
print LINP2 �from modeller import *�;
print LINP2 �from modeller.scripts import complete_pdb�;
print LINP2 ��;
print LINP2 �log.verbose()    \# request verbose output�;
print LINP2 �env = environ()�;
print LINP2 �env.libs.topology.read(file='$(LIB)/top_heav.lib') \# read topology�;
print LINP2 �env.libs.parameters.read(file='$(LIB)/par.lib') \# read parameters�;
print LINP2 ��;
print LINP2 �for i in range(1, 11):�;
print LINP2 �    \# read model file�;
print LINP2 �    code = "Target.BL%04d0001.pdb" % i�;
print LINP2 �    mdl = complete_pdb(env, code)�;
print LINP2 �    s = selection(mdl)�;
print LINP2 �    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='Target.profile',�;
print LINP2 �                  normalize_profile=True, smoothing_window=15)�;
      
#evaluate model
open(LINP3, $modeller_inp8) if $modeller_inp8;
print LINP3 �from modeller import *�;
print LINP3 �from modeller.scripts import complete_pdb�;
print LINP3 ��;
print LINP3 �log.verbose()    # request verbose output�;
print LINP3 �env = environ()�;
print LINP3 �env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology�;
print LINP3 �env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters�;
print LINP3 ��;
print LINP3 �# directories for input atom files�;
print LINP3 �env.io.atom_files_directory = './:../atom_files'�;
print LINP3 ��;
print LINP3 �# read model file�;
print LINP3 �mdl = complete_pdb(env, 'TvLDH.BL00080001.pdb')�;
print LINP3 ��;
print LINP3 �s = selection(mdl)�;
print LINP3 �s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',�;
print LINP3 �              normalize_profile=True, smoothing_window=15)�;

}