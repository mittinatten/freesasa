
use strict;
use warnings;

my @ref_files = `ls reference/*`;
my @sr_files = `ls SR/*/*`;
my @lr_files = `ls LR/*/*`;

my %ref_data;
my %sr_data;
my %lr_data;


read_file(\%ref_data,$_) for @ref_files;
read_file(\%sr_data,$_)  for @sr_files;
read_file(\%lr_data,$_)  for @lr_files;

for my $protein (keys %ref_data) {
    my $r = $ref_data{$protein}[0]->{total};
    my $n = $ref_data{$protein}[0]->{n_atoms};
    for (sort { $a->{par} <=> $b->{par} } @{$sr_data{$protein}}) {
	print "SR $protein $n $r ";
	print "$_->{par} $_->{total} $_->{t}\n";
    }
    for (sort { $a->{par} <=> $b->{par} } @{$lr_data{$protein}}) {
	print "LR $protein $n $r ";
	print "$_->{par} $_->{total} $_->{t}\n";
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

