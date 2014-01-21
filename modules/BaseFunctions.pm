package BaseFunctions;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(stdeviation average);
$VERSION=1.0;

sub stdev{
my ($nn, @values)=@_;
my $sum=0;
for my $i (0..$nn-1){
$sum=$sum+$values[$i];
}
my $avg_val=$sum/$nn;
my $sum_ii=0;
for my $i (0..$nn-1){
my $ii=($values[$i]-$avg_val)^2;
$sum_ii=$sum_ii+$ii;
}
return sqrt($sum_ii/$nn);
}

sub avg{
my ($nn, @values)=@_;
my $sum=0;
for my $i (0..$nn-1){
$sum=$sum+$values[$i];
}
return $sum/$nn;
}
