#!/usr/bin/perl -w
use strict;



my $blast_list = $ARGV[0];
my $kog_dataframe = $ARGV[1];

my %hash_table;

open (LIST, $blast_list) or die "Unable to open input file, $blast_list";
open (KOG, $kog_dataframe) or die "Unable to open input file, $kog_dataframe";

while ( <KOG> )
{
    my $line = $_;
    my @fields = split (/\t/, $line);
    
    if ($fields[0] =~ m/(.+)/)
    {
        $hash_table{$1} = $fields[1];

    }
}

while ( <LIST> )
{
    my $line = $_;
    chomp $line;
    my @fields = split (/\t/, $line);
    my $reads = 0;
    if ($fields[0] =~ m/(.+)/)
        {
            my $blasthit = $1;
            
            if ($hash_table{$blasthit})
            {
                    print "$taxa[$i]", "\t"; 
                }
            }
		
        }
}

close LIST;
close KOG;

