#!/usr/bin/perl

# Meant to be run in this directory to calculate the reference SASA
# values based on the PDB files stored here.

use strict;

my @pdb = `ls *.pdb`;
print "struct residue_sasa sasa_ref[20] = \{\n";
foreach my $p (@pdb) {
    chomp $p;
    my %res;
    my @data = `../src/freesasa $p -n 1000 -R --select="bb, resi 2 and name c+n+o+ca" --select="sc, resi 2 and not name c+n+o+ca" --select="pol2, resi 2 and not symbol c" --select="apol2, resi 2 and symbol c"`;
    foreach (@data) {
        if (/^SEQ A\ +2 (\w\w\w) :\ + (\d+.\d+)/) {
            $res{'name'} = $1;
            $res{'total'} = $2;
        }
        if (/^bb :\ +(\d+.\d+)/) {
            $res{'bb'} = $1;
        }
        if (/^sc :\ +(\d+.\d+)/) {
            $res{'sc'} = $1;
        }
        if (/^pol2 :\ +(\d+.\d+)/) {
            $res{'pol'} = $1;
        }
        if (/^apol2 :\ +(\d+.\d+)/) {
            $res{'apol'} = $1;
        }
    }
    print "    \{.name = \"$res{'name'}\", .total = $res{'total'}, .main_chain = $res{'bb'}, ".
        ".side_chain = $res{'sc'}, .polar = $res{'pol'}, .apolar = $res{'apol'}\},\n";
}
print "\};\n";

