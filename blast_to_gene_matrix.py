#!/usr/bin/env python
#
# blast_to_gene_matrix.py v1 created by WRF 2017-02-10

'''blast_to_gene_matrix.py v1.3 2018-02-26

    convert tabular blast output into matrix of samples to gene counts
    like:
        PM5     PM30      PM70    PM90   PM159
  Gene1  #       #         #       #      #
  Gene2  #       #         #       #      #

    sample information is extracted from the 1st term in the query
  PM159_NoIndex_L004_R1_001_(paired)_trimmed_(paired)_contig_6_[cov=14968]

    sample is PM159

blast_to_gene_matrix.py -i metagenome_blastx_output.tab > sample_gene_matrix.tab

    to add gene description and species
    use -f to include a GenBank style FASTA file

blast_to_gene_matrix.py -i metagenome_blastx_output.tab -f SEED_db.fasta > sample_gene_matrix.tab

    to add signalP hits as .csv file of genename,
blast_to_gene_matrix.py -i metagenome_blastx_output.tab -f SEED_db.fasta -p metagenome_signalp.csv > sample_gene_matrix.tab

'''

import sys
import argparse
import time
import re
from collections import defaultdict

def read_cazy_fasta(cazy_fasta):
	entry_counter = 0
	gi_to_cazy_number = {} # keys are NCBI accessions, values are CAZY GH numbers
	print >> sys.stderr, "# Reading {}".format(cazy_fasta), time.asctime()
	for line in open(cazy_fasta, 'r'):
		if line[0]==">":
			entry_counter += 1
			#gi|474338|gb|AAA17751.1| alpha-amylase [AAA17751.1;GH13
			# gene                    desc           annotation gh_number
			gene, desc = line[1:].strip().split(" ",1) # gene is everything before first space
			annotation, cazy_id = desc.split("[",1)
			gh_number = cazy_id.split(";")[1]
			gi_to_cazy_number[gene] = [gh_number,annotation]
			if entry_counter < 2:
				print >> sys.stderr, "Entries as {}:{} {}".format(gene, gh_number, annotation)
	print >> sys.stderr, "# Counted {} CAZY entries".format(entry_counter), time.asctime()
	return gi_to_cazy_number

def read_cazy_blast(cazy_blast_table):
	hit_counter = 0
	gene_to_gi = {} # key is gene accession ID, value is NCBI GI number blast hit
	print >> sys.stderr, "# Reading {}".format(cazy_blast_table), time.asctime()
	for line in open(cazy_blast_table, 'r'):
		line = line.strip()
		if line and line[0]!="#":
			hit_counter += 1
			lsplits = line.split("\t")
			queryid = lsplits[0]
			subjectid = lsplits[1]
			gene_to_gi[queryid] = subjectid
	print >> sys.stderr, "# Counted {} blast hits against CAZY".format(hit_counter), time.asctime()
	return gene_to_gi

def read_fasta_to_description(fastafile):
	speciesdict = {}
	descdict = {}
	desc_counter = 0
	print >> sys.stderr, "# Reading {}".format(fastafile), time.asctime()
	for line in open(fastafile,'r'):
		if line[0]==">":
			desc_counter += 1
			if desc_counter % 20000 == 0:
				sys.stderr.write(".")
			if desc_counter % 1000000 == 0:
				print >> sys.stderr, "{}".format(desc_counter), time.asctime()
			gene, desc = line[1:].split(" ",1)
			annotation, species = desc.split("[",1)
			species = species.replace("]","").strip()
			speciesdict[gene] = species
			descdict[gene] = annotation.strip() if annotation.strip() else "None"
	print >> sys.stderr, "# Counted {} FASTA descriptions".format(len(descdict)), time.asctime()
	return speciesdict, descdict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', nargs="*", help="tabular blast output")
	parser.add_argument('-d','--delimiter', default="_", help="delimiter for splitting the sample ID [_]")
	parser.add_argument('-s','--similarity', default=60.0, type=float, help="amino acid similarity cutoff [60.0]")
	parser.add_argument('-l','--length', default=50.0, type=float, help="alignment length cutoff [50.0]")
	parser.add_argument('-f','--fasta', help="optional FASTA file to include gene descriptions")
	parser.add_argument('-c','--cazy-blast-results', help="blast results against CAZY enzymes, queries for this table should be subject hits from -i")
	parser.add_argument('-C','--cazy-fasta', help="FASTA format file of CAZY enzymes")
	parser.add_argument('-p','--signalp', help="optional SignalP summarized results file")
	args = parser.parse_args(argv)

	### PARSE OPTIONAL SIGNALP OUTPUT AS tab delimited file OF GENE NAMES
	if args.signalp:
		hassignal = {}
		for line in open(args.signalp,'r'):
			protid = line.split(",")[0]
			hassignal[protid] = "True"

	### READ INPUT DATA
	genecounts = defaultdict( lambda: defaultdict(int) ) # dict of dicts, keys are subject IDs, values are counts
	genelist = {} # store genes as boolean dictionary
	allhits = 0
	hitcounter = 0
	for blastfile in args.input:
		print >> sys.stderr, "# Reading {}".format(blastfile), time.asctime()
		for line in open(blastfile,'r'):
			line = line.strip()
			if line and line[0]!="#":
				allhits += 1
				if float(line.split("\t")[2]) >= args.similarity and float(line.split("\t")[3]) >= args.length:
					hitcounter += 1
					lsplits = line.split("\t")
					queryid = lsplits[0]
					sampleid = queryid.split(args.delimiter)[0]
					try:
						coverage = int( re.search( "\[cov=(\d+)\]", queryid ).group(1) )
					except AttributeError:
						print >> sys.stderr, "NO COVERAGE {}".format(queryid)
						coverage = 1
					subjectid = lsplits[1]
					genelist[subjectid] = True
					genecounts[sampleid][subjectid] += coverage
					if args.signalp:
						if hassignal.get(queryid,False):
							hassignal[subjectid] = "True"
		print >> sys.stderr, "# Counted {} good BLAST hits, of {}".format(hitcounter, allhits), time.asctime()
	print >> sys.stderr, "# Found {} total genes".format(len(genelist))

	### PARSE OPTIONAL FASTA TO GET GENE DESCRIPTIONS

	if args.fasta:
		speciesdict, descdict = read_fasta_to_description(args.fasta)

	### PARSE OPTIONAL CAZY BLAST HITS AND FASTA
	
	gi_to_cazy_number = read_cazy_fasta(args.cazy_fasta) if args.cazy_fasta else {}
	gene_to_gi = read_cazy_blast(args.cazy_blast_results) if args.cazy_blast_results else {}

	### PRINT MATRIX

	# header line
	lineheaders = ["Gene"]
	if args.fasta: # include gene description and species
		lineheaders.extend(["Description","Species"])
	if args.signalp:
		lineheaders.append("Has_signal")
	if args.cazy_blast_results and args.cazy_fasta:
		lineheaders.extend(["CAZY GH","CAZY description"])
	print >> sys.stdout, "\t".join(lineheaders + sorted(genecounts.keys()))

	# print each gene
	for gene in sorted(genelist.keys()):
		countlist = [gene] # build from empty list with only gene name
		if args.fasta:
			countlist.extend([descdict.get(gene,"NONE"), speciesdict.get(gene,"NONE")])
		if args.signalp:
			countlist.append( hassignal.get(gene,"False") )
		if args.cazy_blast_results and args.cazy_fasta:
			ginumber = gene_to_gi.get(gene,None)
			if ginumber is None:
				print >> sys.stderr, "no GI for {}".format(gene)
			countlist.extend( gi_to_cazy_number.get( ginumber , ["NONE","NONE"]) )
		for sk in sorted(genecounts.keys()):
			countlist.append( str(genecounts[sk][gene]) )
		print >> sys.stdout, "\t".join(countlist)
	print >> sys.stderr, "# Process completed", time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
