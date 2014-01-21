#!/usr/bin/perl -w
if ($#ARGV != 2) {die "Program require arguments! [options] [ligand] [receptor fasta]\n";}
my $rootdir=`pwd`;
chomp $rootdir;

use warnings;
use strict;
use lib '/home/wiktor/Zabawki/7TMpipe/TMpipe_mod';
use File::Copy "cp";
#use 7TMpipe_mod::data_library;
use TargetSeqJobs;
use GetReceptorModel;
use LigandStructure;
use RunDocking;
use AnalyzePoses;
#use base_functions;

#read options
open(OPCJE, "< $ARGV[0]") or die "Can not open an input file: $!";
my @param=<OPCJE>;
close (OPCJE);
chomp @param;
my (%params,@para,$lin);
foreach my $lin(@param){
  unless (($lin=~ /^\#/) or ($lin=~ /^$/)) {
    @para=split(/\s+/,$lin);
    $params{$para[0]}=$para[1];
  }
}

my $czas=localtime(time());
my $me=getlogin();
my ($defwyn);	
my $roll=int(rand 10000000) +1;
my $wyniki="run".$roll."out";

#generic output
open(GOUT,"> $wyniki") or die "Can’t write output file: $!";
printf GOUT "Main ouput file generated with perl 7TMpipeline (W.Jurkowski) \n";
printf GOUT "User: %s Time: %s\n", $me, $czas;
printf GOUT "Selected parameters were used: \n";
while ( my ($key, $value) = each(%params) ) {
  printf GOUT "$key = $value\n";
}
printf GOUT "description of output files:\n";

#some settings and inputs
#recpetor and ligand names
my $exe_dir=$params{'bin'};#third party soft installation path
my $ligand=$ARGV[1];#ligand pdb
my $receptor=$ARGV[2];#receptor fasta
my $ind=index($ligand,".");
my $ligand_name=substr($ligand,0,$ind);
$ind=index($receptor,".");
my $receptor_name=substr($receptor,0,$ind);
my $receptor_pir=$receptor_name.".ali";
my $ncore = $receptor_name.".".$ligand_name;

#working directories
my $dir_ligand=$ligand_name."_prep";
my $dir_receptor=$receptor_name."_prep";
my $docking_ad=$receptor_name."_".$ligand_name."_ad";
mkdir($dir_ligand, 0755) if (! -d "$dir_ligand");
mkdir($dir_receptor, 0755) if (! -d "$dir_receptor");
mkdir($docking_ad, 0755) if (! -d "$docking_ad");

#receptor build parameters
my $template_pdbs = $params{'matryce'};
my @templates = split(",",$template_pdbs);
my $ch_pdbs = $params{'chids'};
my @templates_chid= split(",",$ch_pdbs);
my $receptor_model=$params{'receptor_model'};#need to be defined  if stage=4
my $normal_loop=$params{'loop'};
my $dope_loop=$params{'dope_loop'};
my $nr_loop_models=$params{'nr_loop_models'};
my $loop_starts=$params{'loop_starts'};
my @lsegs=split(",",$loop_starts);
my $loop_ends=$params{'loop_ends'};
my @lsege=split(",",$loop_ends);
my $detect_loops=$params{'detect_loops'};#TODO? automatic detection nased on ss prediction
my $ss_starts=$params{'ss_starts'};
my @disulphs=split(",",$ss_starts);
my $ss_ends=$params{'ss_ends'};
my @disulphe=split(",",$ss_ends);
my $scf=$params{'scoring_funct'};
my @core_scorings=();
my ($sel_templ_name,$sel_templ_chid,$n_sel_templ);
my ($tmsegs,$tmsege,$intsegs,$intsege,$extsegs,$extsege);#variables needed for automatic detection of TM helices

#docking parameters
my $ad_grd_npx=$params{'ad_grd_npx'};
my $ad_grd_npy=$params{'ad_grd_npy'};
my $ad_grd_npz=$params{'ad_grd_npz'};
my $ad_grd_space=$params{'ad_grd_space'};
my @bindingsite=$params{'site_center'};
my $dock_flex=$params{'dock_flex'};

#copy initial files
unless($receptor_model eq "undefined"){
	cp($receptor_model,$dir_receptor);
}

#runs stage-wise
if($params{'run_stage'} == 0 or $params{'run_stage'} == 1){
  chdir "$dir_receptor" or die "Can't enter $dir_receptor: $!\n";
  my ($templ_name,$templ_ev,$templ_bit,$templ_ident,$templ_length,$n_templ)=get_templates(\$receptor);#retrieve templates from pdb
  my @templ_prop_ident=@{$templ_ident};
  my @templ_prop_name=@{$templ_name};
  my @templ_pro_ev=@{$templ_ev};
  my @templ_prop_bit=@{$templ_bit};
  my @templ_prop_length= @{$templ_length};
  my $numb_templ=$$n_templ;
  ($sel_templ_name,$sel_templ_chid,$n_sel_templ)=select_template(\@templ_prop_name,\@templ_prop_ident,\@templ_prop_length,\$numb_templ);#selection of top10 single templates
  
  run_topcons(\$receptor,\$exe_dir);#finds TM definitions
  #($tmsegs,$tmsege,$intsegs,$intsege,$extsegs,$extsege)=detect_TM();#parse TM definitions
  #make_pir(\$receptor,\$receptor_name);#saves pir formatted sequence
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 2){
  my ($best,$minim,$cinimini,$cinibest,$cinipdbid);
  chdir "$dir_receptor" or die "Can't enter $dir_receptor: $!\n";
  for my $i (0..$#templates){
    my $tpdbid=$templates[$i];
    my $tchid=$templates_chid[$i];
    make_core_single_inputs(\$tpdbid,\$tchid,\$ncore,\$receptor_pir,\$receptor_name);
    make_eval_inputs(\$ncore,\$nr_loop_models);
    make_model_core(\$ncore);
    #@core_scorings=make_scoring();
    #($best,$minim)=select_core(\@core_scorings,\$scf);
    #`cp "Target.B9999000".$$best.".pdb" $i."T".$i."_".$tpdbid."-".$$best.".pdb"`;	
    #if($cinimini < $$minim){
#	$cinibest=$$best;
#	$cinipdbid=$tpdbid;
#    	$receptor_model=$receptor_name."-T".$i."_".$cinipdbid."-".$cinibest.".pdb";
#    }	
  }
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 3){
  my ($best,$minim);
  chdir "$dir_receptor" or die "Can't enter $dir_receptor: $!\n";
  make_core_inputs(\@templates,\@templates_chid,\$ncore,\$receptor_pir,\$receptor_name);
  make_eval_inputs(\$ncore,\$nr_loop_models);
  make_model_core_m(\$ncore);
  @core_scorings=make_scoring();
  ($best,$minim)=select_core(\@core_scorings,\$scf);
  `cp "Target.B9999000".$$best.".pdb" $receptor_name."-model-".$$best.".pdb"`;	
  $receptor_model=$receptor_name."-model-".$$best.".pdb";
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 4){
  my ($cys1,$cys2);
  my $ss=0;
  chdir "$dir_receptor" or die "Can't enter $dir_receptor: $!\n";
  #for my $i (0..$#{$extsegs}){
  for my $i (0..$#lsegs){
    my $lewy=$lsegs[$i];
    my $prawy=$lsege[$i];
    for my $k(0..$#disulphs){
	if(($disulphs[$k] >=$lewy) and ($disulphe[$k] <=$prawy)){
		$cys1=$disulphs[$k];
		$cys2=$disulphe[$k];
		$ss = 1;
		print "cycy $ss\n";
	}	
    } 
    make_loops_inputs(\$ncore,\$normal_loop,\$dope_loop,\$ss,\$receptor_model,\$lewy,\$prawy,\$nr_loop_models,\$cys1,\$cys2);
    make_eval_inputs(\$ncore,\$nr_loop_models,\$lewy,\$prawy);
    make_model_loops(\$ncore);
    #select_loops(\$ncore,\$normal_loop,\$dope_loop);
  }
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 5){
  chdir "$dir_receptor" or die "Can't enter $dir_receptor: $!\n";
  build_model();
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 6){#ligand
  chdir "$dir_ligand" or die "Can't enter $dir_ligand: $!\n";
  make_ligand(\$ligand_name);
  #optimize_ligand(\$ligand_name);
  #get_partial_charges(\$ligand_name);
  #test_ligand (\$ligand_name);
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 7){#docking
  chdir "$docking_ad" or die "Can't enter $docking_ad: $!\n";
  autdock_inputs(\$ligand,\$receptor,\$ad_grd_npx,\$ad_grd_npy,\$ad_grd_npz,\$ad_grd_space,\@bindingsite,\$dock_flex);
  make_autodock(\$receptor,\$ligand);
  make_poses();
  chdir "$rootdir";
}

if($params{'run_stage'} == 0 or $params{'run_stage'} == 8){#final selection
  chdir "$docking_ad" or die "Can't enter $docking_ad: $!\n";
  chdir "$rootdir";
}

my $koniec=localtime(time());
printf GOUT "Run completed. Time: $koniec\n";
close (GOUT);
