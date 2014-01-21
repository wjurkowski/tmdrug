package LigandStructure;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(make_ligand optimize_ligand_gamess get_partial_charges test_ligand);
$VERSION=1.0;

sub make_ligand{
my ($ligand_name,$opt_method)=@_;
my $ligand=$ligand_name."sdf";
my $ligandmol2=$ligand_name."mol2";
`babel -isdf $ligand -omol2 $ligandmol2`;
}

sub optimize_ligand_gamess{
my ($ligand_name)=@_;
my $gamess_in=$ligand_name."in";
my $ligand=$ligand_name."sdf";
open(HEAD, "> header") or die "Can not open an input file: $!";
print HEAD "'$CONTRL  ICHARG=0 MPLEVL=0'\n";
print HEAD "         RUNTYP=OPTIMIZE SCFTYP=RHF EXETYP=RUN\n";
print HEAD "          MULT=1 UNITS=ANGS MAXIT=200\n";
print HEAD "          COORD=CART                           $END\n";
print HEAD "$SCF     CONV=1.0E-08                         $END\n";
print HEAD "$SYSTEM  TIMLIM=50000 MWORDS=10 MEMDDI=0      $END\n";
print HEAD "$STATPT  NSTEP=200 OPTTOL=1.0E-06\n";
print HEAD "         HESS=CALC IHREP=10                   $END\n";
print HEAD "$FORCE   METHOD=ANALYTIC VIBANL=.F.           $END\n";
print HEAD "$BASIS   GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0\n";
print HEAD "         DIFFSP=.F.                           $END\n";
print HEAD "$GUESS   GUESS=HUCKEL                         $END\n";

`babel -isdf $ligand -ogukin $gamess_in`;
`cat header $gamess_in > g_opt.in`;
`rungamess g_opt.in g_opt.out`;
`rungamess g_charges.in g_charges.out`;
}

sub get_antechamber_charges{
my ($ligand_name,$charge_type)=@_;
my $ligand=$$ligand_name."pdb";
my $ligand_out=$$ligand_name."mol2";
`antechamber -i $ligand -fi pdb -o $ligand_out -fo mol2 -c bcc -s 2`; 
}

sub run_RED{
my ($ligand_name,$gamess_run)=@_;
my $ligand=$$ligand_name."pdb";
`perl Ante_RED.pl $ligand > Ante_RED.log`;
#add "REMARK TITLE Dimethylalanine-dipeptide" jezeli nie ma w pliku p2n
#check charge and multiplicity
#check connectivity
#check atom names
if ($$gamess_run == 0){
`perl RED-vIII-on.pl > RED-vIII.log`; 
}
elsif ($$gamess_run == 1){
optimize_ligand_gamess($ligand_name);
`perl RED-vIII-off.pl > RED-vIII.log`; 
}
}

sub test_ligand{
my ($ligand_name)=@_;
print "pustka\n";
}



