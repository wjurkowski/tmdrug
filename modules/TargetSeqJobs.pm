package TargetSeqJobs;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::AlignIO;
use Bio::SearchIO;
use File::Copy;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_templates select_template detect_TM run_topcons make_pir);
$VERSION=1.0;

sub get_templates{
my ($receptor)=@_;
my (@templ_name,@templ_ev,@templ_bit,@templ_ident,@templ_length);
my $rec=$$receptor;
my @res=`blastp -task blastp -query $rec -db pdbaa -outfmt '6 sacc evalue bitscore pident length'`;
$n=0;
foreach my $lin(@res){
  my @linia=split(/\s+/,$lin);
  if($linia[1] <= 1e-05){  
 	$n++;
 	$templ_name[$n]=$linia[0];
 	$templ_ev[$n]=$linia[1];
 	$templ_bit[$n]=$linia[2];
 	$templ_ident[$n]=$linia[3];
	$templ_length[$n]=$linia[4];
  }
}
if ($n > 0){
 return (\@templ_name,\@templ_ev,\@templ_bit,\@templ_ident,\@templ_length,\$n);
}
}

sub select_template{
my ($templ_name,$templ_ident,$templ_length,$n_templ)=@_;
my (@best,@second,@sel_templ,@sel_templ2,@sel_templ_name,@sel_templ_chid);
$best[0]=0;
my $m=0;
my $k=0;
for my $i(2..$$n_templ){
  if($$templ_ident[$i] >= $$templ_ident[1]-10){
    $m++;
    $best[$m]=$i;
  }
  else{
    $second[$k]=$i;
    $k++;
  }
}
my $najdl=0;
for my $i(1..$#best){
#print "leheleh $i $$templ_ident[$i]\n";
  if($$templ_length[$best[$i]] >$najdl+$najdl*0.1){
    $najdl=$best[$i];
    unshift (@sel_templ,$best[$i]);
  }
}

$najdl=0;
for my $i(1..$#second){
  if($$templ_length[$second[$i]] >$najdl+$najdl*0.1){
    $najdl=$second[$i];
    unshift (@sel_templ2,$second[$i]);
  }
}
push (@sel_templ,@sel_templ2);
my $n=10;
if ($#sel_templ < 10) {$n=$#sel_templ;}
  for my $i (0..$n) {
    ($sel_templ_name[$i],$sel_templ_chid[$i])=split(/-/,$$templ_name[$sel_templ[$i]]);
  }
  return (\@sel_templ_name,\@sel_templ_chid,\$n);
}

sub detect_TM{
my (@tmsegs,@tmsege,@intsegs,@intsege,@extsegs,@extsege);
open(TOPC, "<$topcons_out") or die "Cant open output file: $! $topcons_out"; #open top cons output
my @file=<TOPC>;
close (TOPC);
chomp @file;

my $n=1;
my @topol=();
DUPA: for my $i (0..$#file) {
  if ($file[$i]=~/^"TOPCONS predicted"/){
	for my $j ($i+1..$#file) {
	  last DUPA if $file[$j]=~/^\s/;
	  my @linia=split(//,$file[$j]);
	  push(@topol,$linia);  
	}
  }
}

$tm=0;
$in=0;
$ex=0;
if ($topol[0] eq "M"){
  $tm++;
  $tmsegs[$tm]=1;
}
if ($topol[0] eq "i"){
  $in++;
  $intsegs[$in]=1;
}
if ($topol[0] eq "o"){
  $ex++;
  $extsegs[$ex]=1;
}

for $aa (1..$#topol){
  if ($topol[$aa] eq "M" and $topol[$aa-1] eq "i"){
	$tm++;
	$tmsegs[$tm]=$aa; 
	$intsege[$in]=$aa;
  }
  if ($topol[$aa] eq "M" and $topol[$aa-1] eq "o"){ 
	$tm++;
	$tmsegs[$tm]=$aa; 
	$extsege[$ex]=$aa;
  }
  if ($topol[$aa] eq "i" and $topol[$aa-1] eq "M"){ 
  $tmsege[$tm]=$aa; 
  $in++;
  $intsegs[$in]=$aa;
  }
  if ($topol[$aa] eq "o" and $topol[$aa-1] eq "M"){ 
  $tmsege[$tm]=$aa; 
  $ex++;
  $extsegs[$ex]=$aa;
  }
}

return (\@tmsegs, \@tmsege, \@intsegs, \@intsege, \@extsegs, \@extsege);
}

#topcons
sub run_topcons{
print "\n";
print "TOPCONS\n";
my ($receptor,$exe_dir)=@_;
my $wyjscie="topcons_works";
my $rec=$$receptor;
mkdir($wyjscie, 0755) if (! -d "$wyjscie");
copy("$rec","$wyjscie/$rec");
chdir "$wyjscie" or die "Can't enter $wyjscie: $!\n";
my $target=substr($rec,0,index($rec,"."));
my $recl=$target.".fa";
my $protnamef="protein_list";
open(PRNF, ">$protnamef") or die "Can’t open output file: $! protein_list";
print PRNF "$target\n";
close(PRNF);

run_spoctopus(\$rec,\$exedir);
run_prodiv(\$rec,\$exedir);
run_scampi(\$rec,\$exedir);
#`TOPCONS.sh $protnamef $wyjscie $spoctopusdir $prodivdir $scampidir`;
chdir "../";
}

#spoctopus
sub run_spoctopus{
print "\n";
print "SPOCTOPUS\n";
my ($receptor,$exe_dir)=@_;
my $target=substr($rec,0,index($$receptor,"."));
my $recl=$target.".fa";
my $protnamef="protein_list";
open(PRNF, ">$protnamef") or die "Can’t open output file: $! protein_list";
print PRNF "$target\n";
close(PRNF);

my $spoctopusdir="_spoctopus";
my $pssm_prf_dir="pssm_prf_dir";
my $raw_prf_dir="raw_prf_dir"; 
my $prf_work_dir="prf_work_dir"; 
my $prf_work_out_dir="prf_work_out_dir"; 
my $modhmm=$$exe_dir."/spoctopus/modhmmblast";
mkdir("$spoctopusdir", 0755) if (! -d "$spoctopusdir");
mkdir("$spoctopusdir/$pssm_prf_dir", 0755) if (! -d "$spoctopusdir/$pssm_prf_dir");
mkdir("$spoctopusdir/$raw_prf_dir", 0755) if (! -d "$spoctopusdir/$raw_prf_dir");
mkdir("$spoctopusdir/$prf_work_dir", 0755) if (! -d "$spoctopusdir/$prf_work_dir");
mkdir("$spoctopusdir/$prf_work_out_dir", 0755) if (! -d "$spoctopusdir/$prf_work_out_dir");

chdir "$spoctopusdir" or die "Can't enter $spoctopusdir: $!\n";
#`blastp -task blastp -outfmt 6 -evalue 1.e-5 -query $rec -db pdbaa -num_alignments 250 -out $raw_prf_dir/$target."msa"`;
copy("../$rec","$recl");
copy("../$protnamef","$protnamef");
`blastpgp -m 6 -F F -e 1.e-5 -i $recl -d pdbaa -o $raw_prf_dir/$target.msa`;
#run msa2prf.sh: 
`perl $modhmm/msa62mod.pl $protnamef 0 $raw_prf_dir $prf_work_dir .`;
`perl $modhmm/mod2modquery.pl $protnamef $prf_work_dir $prf_work_out_dir`;
`perl $modhmm/mod_upper_caser.pl $protnamef $prf_work_out_dir`;
`perl $modhmm/mod2prfmod_nolabelfile.pl $protnamef . $prf_work_out_dir $raw_prf_dir`;
`rm -rf $prf_work_dir`;
`rm -rf $prf_work_out_dir`;
`cp $recl $pssm_prf_dir/$recl`;
`cp $recl $pssm_prf_dir/$target."chd"`;
chdir "$pssm_prf_dir";
`blastpgp -j 2 -m 6 -F F -e 1.e-5 -i $recl -d pdbaa -C $target."chk" >psiblast.out`;
#`psiblast -num_iterations 2 -outfmt 6 -evalue 1.e-5 -query $$receptor -db pdbaa -out_ascii_pssm $pssm_prf_dir/$target.".chk" -out psiblast.out`;
`ls -1 *.chk >DATABASE.pn`;
`ls -1 *.chd >DATABASE.sn`;
`makemat -S 1 -P DATABASE`;
`perl $modhmm/mtx2prf.pl $target."mtx" .`;
chdir "../";
my $spocti=$$exe_dir."/spoctopus/SPOCTOPUS.sh";
`bash $spocti $protnamef $pssm_prf_dir $raw_prf_dir .`;
chdir "../";
}

#prodiv
sub run_prodiv{
print "\n";
print "PRODIV\n";
my ($receptor,$exe_dir)=@_;
my $target=substr($rec,0,index($$receptor,"."));
my $recl=$target.".fa";
my $protnamef="protein_list";
open(PRNF, ">$protnamef") or die "Can’t open output file: $! protein_list";
print PRNF "$target\n";
close(PRNF);

my $prodivdir="_prodiv";
my $prodiv_feed="prodiv_in";
mkdir($prodivdir, 0755) if (! -d "$prodivdir");
mkdir("$prodivdir/$prodiv_feed", 0755) if (! -d "$prodivdir/$prodiv_feed");
copy("$rec","$prodivdir/$recl");
chdir "$prodivdir" or die "Can't enter $prodivdir: $!\n";
#`blastp -task blastp -evalue 1.e-5 -query $rec -db swissprot -out $prodiv_feed/$rec.".blast"`;
`blastpgp -e 1.e-5 -i $recl -d nr -o $prodiv_feed/$target.blast`;
chdir "$prodiv_feed" or die "Can't enter $prodiv_feed: $!\n";
my $blast_report = new Bio::SearchIO ('-format' => 'blast', '-file'   => $target.".blast");
print "dupa $blast_report\n";
my $result = $blast_report->next_result;
my $out = Bio::AlignIO->newFh(-format => 'selex' );
open(PRNF, ">$b_rep") or die "Can’t open output file: $! protein_list";
while( my $hit = $result->next_hit()) {
	while( my $hsp = $hit->next_hsp()) { 
	    my $aln = $hsp->get_aln();
#	    print $out $aln;
	}
}
print "dupa $rec\n";
chdir "../";
#ewentualnie najpierw zapisac faste 
#my $in  = Bio::AlignIO->newFh(-file => $inputfilename ,'-format' => 'fasta');
#my $out = Bio::AlignIO->newFh('-format' => 'selex');
#print $out $_ while <$in>;
`selex2mod.pl $protnamef $prodiv_feed $prodivdir`;
`all_tmhmm_runner.pl prodiv $prodiv_feed $prodivdir`;
chdir "../";
}

#scampi
sub run_scampi{
print "\n";
print "SCAMPI\n";
my ($receptor,$exe_dir)=@_;
my $target=substr($rec,0,index($$receptor,"."));
my $recl=$target.".fa";
my $protnamef="protein_list";
open(PRNF, ">$protnamef") or die "Can’t open output file: $! protein_list";
print PRNF "$target\n";
close(PRNF);

my $scampidir="_scampi";
mkdir($scampidir, 0755) if (! -d "$scampidir");
copy($rec,$scampidir/$rec) or die "Copy failed: $!";
chdir "$scampidir" or die "Can't enter $scampidir: $!\n";
`SCAMPI_run.pl $rec 1`;
chdir "../";
}


sub make_pir{
my ($receptor,$receptor_name)=@_;
my ($k1,$k2,$j,$line,@seqline,@aasekw);
my $receptor_pir=$$receptor_name.".ali";
open(SEQ, "<$$receptor") or die "Cant open output file: $! $$receptor";
open(SEQPIR, ">$$receptor_pir") or die "Can’t open output file: $! $receptor_pir";
#read fasta
my @fasta=<SEQ>;
close (SEQ);
chomp @fasta;
foreach $line (@fasta){
  unless($line=~/^>/){@seqline=split(//,$line);}
  push (@aasekw, @seqline);
}
printf SEQPIR ">P1;$$receptor_name";
printf SEQPIR "sequence:$$receptor_name:::::::0.00: 0.00";

my $llin=int $#aasekw/80.0;

for $k(1..$llin+1){# print in lines not longer than 80 characters
	$k1=($k-1)*80+1;
	$k2=$k*80;
	for $j($k1..$k2){
	  if($j <= $#aasekw){
	  printf SEQPIR ("%s",$aasekw[$j]);}
	}
	print SEQPIR "\n";	
}
close (SEQPIR);	
}
1;
