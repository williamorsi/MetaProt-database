#!/usr/bin/perl -w
use strict;

#to run this script you need to files, a list of accessions that had hits to atribacteria (input 2) and a diamond output file (input 1)
my $input_file1 = $ARGV[0];

#this is the COG_ids.txt dataframe file
my $input_file2 = $ARGV[1];

#this is the output from the part1.pl script, with diamond_db_accessions / COG accession / COD id
my %hash_table;


open (INPUT1, $input_file1) or die "Unable to open input file, $input_file1";
open (INPUT2, $input_file2) or die "Unable to open input file, $input_file2";

while ( <INPUT1> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    
    if ($fields[0] =~ m/(\S+)/)
    {
        $hash_table{$1} = $fields[1];
    }
}

while ( <INPUT2> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    if ($fields[2] =~ m/(\S+)/)
        {
            my $kog_accession = $1;
            my $genbank_accession = $fields[0];
            my $kog_description = $fields[2];
             print "$genbank_accession\t", "$kog_accession\t" ;

            if ($hash_table{$kog_accession})
            {
                         
                  print "$hash_table{$kog_accession}", "\n"; 
                    #this will print the genback accession, kog accession, and the KOG description
                }
			
            }

        }
close INPUT1;
close INPUT2;

