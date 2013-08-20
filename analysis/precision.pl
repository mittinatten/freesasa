
use strict;
use warnings;
use Statistics::Descriptive;

my @ref_files = `ls reference/*`;
my @sr_files = `ls SR/*/*`;
my @lr_files = `ls LR/*/*`;

my %ref_data;
my %sr_data;
my %lr_data;
#my %sum;
my %val;

read_file(\%ref_data,$_) for @ref_files;
read_file(\%sr_data,$_)  for @sr_files;
read_file(\%lr_data,$_)  for @lr_files;

for my $protein (keys %ref_data) {
    get_val("SR",$protein,%sr_data);
    get_val("LR",$protein,%lr_data);
}


for (keys %val) {
    #print (join "\n", @{$val{$_}{t}});
    print "STAT $_ ";
    print_quantile (@{$val{$_}{t}});
    print_quantile (@{$val{$_}{d}});
    print "\n";
}

sub print_quantile {
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@_);
    
    # this avoids a strange bug in percentile(), don't understand it fully yet
    my ($x1, $idx1) = $stat->percentile(5); 
    
    printf("%6f %6f %6f %6f %6f ", 
	   $stat->quantile(2), $stat->quantile(1), 
	   $stat->quantile(3), ($stat->get_data)[$idx1],
	   $stat->percentile(95));
    $stat->clear();
}


sub get_val {
    my ($name,$protein,%data) = @_;
    my $r = $ref_data{$protein}[0]->{total};
    my $n = $ref_data{$protein}[0]->{n_atoms};
    for (@{$data{$protein}}) {
	#print "$name $protein $n $r ";
	#print "$_->{par} $_->{total} $t $d\n";
	my $key = "$name\_$_->{par}";
	my $t = $_->{t};
	my $d = abs($r-$_->{total})/$n;
	push @{$val{$key}{t}}, $t/$n;
	push @{$val{$key}{d}}, $d;
    }
}

sub read_file 
{
    my $data = shift;
    my $filename = shift;
    chomp $filename;
    open my $fh, "<", $filename 
	or die "Failed to open file '$filename' for reading: $!";
    my ($pdbfile, $par);
    my %tmp;
    while (<$fh>) {
	next if (/^#/);
	$pdbfile = $1      if (/^File: (\S+)/);
	# both should not be present in the same file
	$tmp{par} = $1     if (/^N_testpoint: (\d+)/);
	$tmp{par} = $1     if (/^d_slice: (\S+)/);
	$tmp{t} = $1       if (/^time_elapsed: (\S+)/);
	$tmp{n_atoms} = $1 if (/^n_atoms: (\d+)/);
	$tmp{total} = $1   if (/^total\s+(\S+)/);
    }
    $pdbfile =~ s/.*\/(\w\w\w\w).pdb/$1/;
    push @{$data->{$pdbfile}}, \%tmp;

}

