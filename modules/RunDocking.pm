package RunDocking;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(autodock_inputs make_autodock make_poses);
$VERSION=1.0;

sub autodock_inputs {
my ($ligand,$receptor,$ad_grd_npx,$ad_grd_npy,$ad_grd_npz,$ad_grid_space,@bindingsite,$dock_flex)=@_;
`prepare_ligand4.py –l $ligand -d "ligand_dictionary.py"`;
`examine_ligand_dict.py > summary.txt`;
if ($dock_flex == 0) {`prepare_receptor4.py -r $receptor -A 'hydrogens'; #require rceptor file name`;}
if ($dock_flex == 1) {`prepare_flexreceptor4.py -r $receptor' -s $bindingsite ; #requires list of residues to be flexible`;}
my $nres=$#bindingsite;
my $sumx=0;
my $sumy=0;
my $sumz=0;

for $i (0..$nres){
$sumx=$sumx+$ikses[$i];
$sumx=$sumy+$uajs[$i];
$sumx=$sumz+$zetas[$i];
}
my $avx=$sumx/$nres;
my $avy=$sumy/$nres;
my $avz=$sumz/$nres;
my $ligandpdbqt=$ligand."qt";
my $receptorqt=$receptor."qt";

`prepare_dpf4.py -l $ligandpdbqt -r $receptorqt`;
`prepare_gpf4.py -l $ligandpdbqt -r $receptorqt -p gridcenter="$avx $avy $avz" -p spacing=$ad_grid_space -p npts=$ad_grd_npx,$ad_grd_npy,$ad_grd_npz`;
}
  
sub make_autodock {
my ($receptor,$ligand)=@_;
my $dpf=$receptor.$ligand."dpf";
my $gpf=$receptor.$ligand."gpf";
my $dlg=$receptor.$ligand."dlg";
my $glg=$receptor.$ligand."glg";
`autogrid4 -p $gpf -l $glg`;
`autodock4 -p $dpf -l $dlg`;
}

sub make_poses {
my ($rmsdtol,$receptor,$ligand,$nposes)=@_;
#`summarize_results4.py –d $d –t $rmsdtol –B –a –o ../etc/summary_2.0.txt`;
my $dlg=$receptor.$ligand."dlg";
`getdocked $dlg`;
}




