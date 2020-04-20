#!/usr/bin/env python
#
# v3.0 switched to reading queries as dictionary 2016-02-17
# v2.5 added file format for fastq sequences 2016-02-10
# v2.4 added fasta mode to find fasta seqs in other dbs 2014-11-25
# v2.3 fixed error when list has empty lines 2014-10-16
# v2.2 added no-find option to ignore string searches 2014-03-27
# v2.1 clarified most exceptions to KeyError 2013/11/14
# v2.0 converted to argparse, to allow complex arguments 2013/08/17
# v1.6 reads in namelist only one time, added quiet mode 2013/08/16
# v1.5 will adapt if second file or list have redundant names 2013/07/24
# v1.4 allows use of multiple files 2013/04/04
# v1.3 corrected find list with ">" 26/03/2013
# v1.2 renamed to more conventional variables, also accept fasta names 20/03/2012
# v1.1 30/03/2012

"""getAinB.py v3.1  last modified 2016-08-24
  takes either a searchterm (such as a contig name, unique tag, or gene name)
  or a file of searchterms (such as a list of contigs, or transcripts)
  and prints all those in a fasta file that match any term

  terms are case-sensitive, even in search mode (-s)

    usage examples:
getAinB.py namelist.txt file1.fasta file2.fasta > found_seqs.fasta
getAinB.py namelist.txt *.fasta > found_seqs.fasta
getAinB.py contig123 contigs.fasta > contig123.fasta

getAinB.py -s DDC_HUMAN uniprot_sprot.fasta > human_ddc.fasta
getAinB.py -s DDC uniprot_sprot.fasta > decarboxylases.fasta

    to print names which are not found, use -v (verbose)
    to only print seqs and no stderr, use -q (quiet)
"""

import os
import sys
import argparse
import time
from Bio import SeqIO

def make_query_dict(query_name, isfasta):
	querydict = {}
	nospacedict = {}
	if os.path.isfile(query_name):
		for line in open(query_name,'r').readlines():
			line = line.strip()
			if not line: # skip empty lines
				continue
			if isfasta and line[0]!=">": # if using fasta file and line is not a header, then skip
				continue
			if line[0]==">": # trim any fasta headers
				line = line[1:]
			# querydict serves as a counter for each key, for found/not found
			querydict[line] = 0
			# nospacedict values are just links to the querydict
			nospacedict[line.split(" ")[0]] = line # the key is the raw term before any spaces
	else: # for single search terms
		querydict[query_name] = 0
		nospacedict[query_name.split(" ")[0]] = query_name
	print >> sys.stderr, "# Searching for {} entries".format( len(querydict) )
	return querydict, nospacedict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('query_name', help="file of names or search term")
	parser.add_argument('search_file', nargs="*", default = '-', help="fasta format files")
	parser.add_argument('--format', default="fasta", help="file format, fasta or fastq")
	parser.add_argument('-d','--delimiter', help="delimiter for splitting sequence IDs")
	parser.add_argument('-f','--fasta', action="store_true", help="assume query is normal fasta file")
	parser.add_argument('-F','--file-name', action="store_true", help="append file name to sequence ID if searching for the same key in many files")
	parser.add_argument('-q','--quiet', action="store_true", help="quiet, will not print stderr messages")
	parser.add_argument('-s','--search', action="store_true", help="search for sequences that do not match exactly, including partial matches, so g1 returns g1, g10, etc.")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose, print names which are not found")
	args = parser.parse_args(argv)

	querydict, nospacedict = make_query_dict( args.query_name, args.fasta )

	for fastafile in [n for n in args.search_file if os.path.isfile(n)]:
		# create a new copy of querydict for each file, so counts are not cumulative
		qdcopy = dict(querydict)
		seqcounter = 0
		if not args.quiet:
			print >> sys.stderr, "# Searching for sequences from {}".format(fastafile), time.asctime()
		for seqrec in SeqIO.parse(fastafile, args.format):
			seqcounter += 1
			# get searchid from seqrec.id and change if needed
			searchid = str(seqrec.id)
			# if searching multiple files, optionally append the filename to the ID
			if args.file_name:
				seqrec.id = "{}__{}".format(seqrec.id, fastafile)
			if args.delimiter: # if a standard delimiter is given, split the ID at that point one time
				searchid = searchid.split(args.delimiter,1)[0]  # thus g123.t1 becomes g123
			# begin search
			if qdcopy.get(searchid, None) is not None:
				qdcopy[searchid] = qdcopy.get(searchid) + 1
				wayout.write(seqrec.format( args.format ))
			elif nospacedict.get(searchid, None): # if not in qdcopy, check for ID without spaces
				qdcopy[nospacedict.get(searchid)] += 1 # value in nospacedict is a key for qdcopy
				wayout.write(seqrec.format( args.format ))
			elif args.search:
				writeseq = False
				for qkey in qdcopy.iterkeys(): # first try string.find() in keys of qdcopy
					if searchid.find(qkey) != -1:
						qdcopy[qkey] = qdcopy.get(qkey) + 1
						writeseq = True
					# do not break so that multiple keys can be counted
				if writeseq: # if any keys were found, then write the sequence
					wayout.write(seqrec.format( args.format ))
				else: # otherwise try string.find() in keys of nospacedict
					for nkey in nospacedict.iterkeys():
						if searchid.find(nkey) != -1:
							qdcopy[nospacedict.get(nkey)] += 1
							writeseq = True
					if writeseq: # as above
						wayout.write(seqrec.format( args.format ))
		found = len([v for v in qdcopy.values() if v > 0])
		multicounts = len([v for v in qdcopy.values() if v > 1])
		notfound = len([v for v in qdcopy.values() if v==0])
		if not args.quiet:
			print >> sys.stderr, "# {} had {} sequences".format(fastafile, seqcounter), time.asctime()
			print >> sys.stderr, "# Found {} names from {} in {}".format(found, args.query_name, fastafile)
			if multicounts:
				print >> sys.stderr, "# {} names from A were found multiple times".format(multicounts)
				if args.verbose:
					for k,v in qdcopy.iteritems():
						if v > 1: # only multicopy queries
							print >> sys.stderr, "{} was found {} times".format(k,v)
			print >> sys.stderr, "# Did not find {} names".format(notfound)
			if args.verbose:
				for k,v in qdcopy.iteritems():
					if v == 0: # only queries not found at all
						print >> sys.stderr, "{} was not found".format(k)
	if not args.quiet:
		print >> sys.stderr, "# Finished", time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
