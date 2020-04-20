#!/usr/bin/perl -w
use strict;



my $COG_mapping = $ARGV[0];
my $matrix = $ARGV[1];

my %hash_table;

open (COGS, $COG_mapping) or die "Unable to open contig/reads input file, $COG_mapping";
open (MATRIX, $matrix) or die "Unable to open matrix input file, $matrix";


while ( <COGS> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    
    if ($fields[0] =~ m/(\S+)/)
    {
        $hash_table{$1} = $fields[2];
    }
}

while ( <MATRIX> )
{
    my $line = $_;
    chomp $line;
    my @fields = split (/\t/, $line);
    my $reads = 0;
    if ($fields[0] =~ m/(\S+)/)
        {
            my $genbank_accession = $1;
            
            
            for (my $i = 0; $i <= $#fields; $i++)
            {
                print "$fields[$i]", "\t";
            }
            
            
            if ($hash_table{$genbank_accession})
            {
                    print "$hash_table{$genbank_accession}"; 
            
                }
                
                else {
                    print "nohit\n"; 
            
                }
            }
            
           
        }



close COGS;
close MATRIX;

