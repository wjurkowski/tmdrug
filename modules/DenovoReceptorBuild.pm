package DenovoReceptorBuild;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::AlignIO;
use Bio::AlignIO;
use File::Copy;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_templates detect_TM run_topcons make_pir);
$VERSION=1.0;


sub build_nonh{
my ($tmsegs, $tmsege)=@_;
print GOUT "build nonhomologous structure\n";
}
