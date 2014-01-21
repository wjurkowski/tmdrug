package AnalyzePoses;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(pose_cluster);
$VERSION=1.0;

sub pose_cluster {
my (@posermsd,@cluster);
my $is_sth_lonely=1;
my $nclus=0;
for my $i (1..$nposes){push(@posermsd,"1");}

foreach my $i (0..$nposes){
	$klastry[$nclus]=@cluster;
	my @noncluster=();
	my @cluster=();
	
	if($is_sth_lonely==1){
	for my $k (0..$#posermsd){
		if ($posermsd[$k] < 0.5) {push (@cluster,$poseid[$k]);}
		else {push (@noncluster,$poseid[$k]);}
	}
	if ($#cluster > 0){$nclus++;}
	if ($#noncluster != 0){
		$is_sth_lonely=1;
		for my $j (0..$#noncluster){`pr_alchem prm-RMSD $pose.$noncluster[0] $pose.$j >>pairs_RMSD`;}
		open(RMS1, "< pairs_RMSD") or die "Can not open an input file: $!";
		my @rmsds=<RMS1>;
		close (RMS1);
		chomp @rmsds;
		my $n=0;
		foreach my $line(@rmsds){
		  my @values=split(/\s+/,$line);
		  $posermsd[$n]=$values[4];
		  $poseid[$n]=$i;
		  $n++;
		}
	}
	else {$is_sth_lonely=0;}
	}
}
return (@klastry);
}
