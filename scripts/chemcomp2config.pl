# Simple script to read an entry from the Chemical Component
# Dictionary (ftp://ftp.wwpdb.org/pub/pdb/data/monomers) and create
# lines for a classifier using the ProtOr classes in the paper by Tsai
# et al. (http://www.ncbi.nlm.nih.gov/pubmed/10388571).

## Example ##
# Using the entry for ALA
#
#  RESIDUE   ALA     13
#  CONECT      N      3 CA   H    H2  
#  CONECT      CA     4 N    C    CB   HA  
#  CONECT      C      3 CA   O    OXT 
#  CONECT      O      1 C   
#  CONECT      CB     4 CA   HB1  HB2  HB3 
#  CONECT      OXT    2 C    HXT 
#  CONECT      H      1 N   
#  CONECT      H2     1 N   
#  CONECT      HA     1 CA  
#  CONECT      HB1    1 CB  
#  CONECT      HB2    1 CB  
#  CONECT      HB3    1 CB  
#  CONECT      HXT    1 OXT 
#  END   
#  HET    ALA             13
#  HETNAM     ALA ALANINE
#  FORMUL      ALA    C3 H7 N1 O2
#
# We get the output
#
#  ALA N N3H2
#  ALA CA C4H1
#  ALA C C3H0
#  ALA O O1H0
#  ALA CB C4H3
#  ALA OXT O2H1
#
# These lines can be added to the atoms-section of a classifier
# config-file, and then the types N3H2, C4H1 need to be defined in the
# types-section. See share/protor.config for a template including
# all the standard amino acids, nucleic acids, plus some common
# non-standard entries. In the setting of a polypeptide chain the
# N-atom should be N3H1, since it is linked to the next residue, but
# since N3H1 and N3H2 have the same radius this doesn't necessarily
# have to be fixed.
#
# The script does not do anything clever, so if your entry contains
# some unusual elements the output should be checked (the atom SE in
# SEC becomes S2H1 for instance).
#

use strict;
my $res;
while (<>) {
    if (/RESIDUE\s+(\w+)/) { $res = $1; next }
    if (/CONECT/) {
        my @fields = split /\s+/, $_;
        shift @fields;
        my $atom = shift @fields;
        next if ($atom =~ m/^H/);
        my $val = shift @fields;
        my $nh = 0;
        foreach (@fields) {
            ++$nh if ($_ =~ m/^H/);
        }
        print "$res $atom ",substr($atom,0,1),$val,"H$nh\n";
    }
}
