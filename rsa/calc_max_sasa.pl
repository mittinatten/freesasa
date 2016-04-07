#!/usr/bin/perl

# Meant to be run in this directory to calculate the reference SASA
# values based on the PDB files stored here. Defines all C atoms as apolar
# and all else as polar, i.e. doesn't work for OONS.

use strict;

my @pdb = `ls *.pdb`;
my @res_array;
print ">>>> C code";
print "freesasa_residue_sasa sasa_ref[20] = \{\n";
foreach my $p (@pdb) {
    chomp $p;
    my %res;
    my @data = `../src/freesasa $p -n 1000 -R --select="bb, resi 2 and name c+n+o+ca" --select="sc, resi 2 and not name c+n+o+ca" --select="pol2, resi 2 and not symbol c" --select="apol2, resi 2 and symbol c"`;
    foreach (@data) {
        if (/^SEQ A\ +2 (\w\w\w) :\ + (\d+.\d+)/) {
            $res{name} = $1;
            $res{total} = $2;
        }
        if (/^bb :\ +(\d+.\d+)/) {
            $res{bb} = $1;
        }
        if (/^sc :\ +(\d+.\d+)/) {
            $res{sc} = $1;
        }
        if (/^pol2 :\ +(\d+.\d+)/) {
            $res{pol} = $1;
        }
        if (/^apol2 :\ +(\d+.\d+)/) {
            $res{apol} = $1;
        }
    }
    push @res_array, \%res;
    print "    \{.name = \"$res{name}\", .total = $res{total}, .main_chain = $res{bb}, ".
        ".side_chain = $res{sc}, .polar = $res{pol}, .apolar = $res{apol}\},\n";
}
print "\};\n";

# The following could potentially be an extension to the config-files
# to let users specify their own RSA configuration

#print "\n>>>> Config-file\n";
#print "reference-sasa:\n";
#foreach my $res (@res_array) {
#    my $name = $res->{name};
#    print "$name total $res->{total}\n";
#    print "$name apolar $res->{apol}\n";
#    print "$name polar $res->{pol}\n";
#    print "$name bb $res->{bb}\n";
#    print "$name sc $res->{sc}\n";
#    print "\n";
#}
