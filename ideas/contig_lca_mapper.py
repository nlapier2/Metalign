import argparse, ast, math, os, subprocess, sys, time


start = time.time()  # start a program timer
__location__ = os.path.realpath(os.path.join(os.getcwd(),
								os.path.dirname(__file__))) + '/'
RANKS = ['superkingdom', 'phylum', 'class', 'order',
		'family', 'genus', 'species', 'strain']


def echo(msg, verbose):
	if not verbose:
		return
	global start
	seconds = time.time() - start
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	hms = "%02d:%02d:%02d" % (h, m, s)
	print('['+hms+'] ' + msg)


def parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Compute abundance estimations for species in a sample.')
	parser.add_argument('infile',
		help='sam file mapping contigs to accession unique regions. Required.')
	parser.add_argument('dbinfo', help = 'Location of db_info file.')
	parser.add_argument('--output', default='abundances.tsv', help='Output file name.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--verbose', action='store_true',
		help='Print verbose output.')
	args = parser.parse_args()
	return args


# Given a taxonomic lineage, get the rank of the taxid
def get_taxid_rank(taxlin):
	# get number of empty tax ranks at end of taxlin
	end_empty, splits = 0, taxlin.split('|')
	for i in range(1, len(taxlin)+1):
		if splits[-i] == '':
			end_empty += 1
			continue
		break
	return RANKS[-(end_empty+1)]


# Parses information in dbinfo file, which maps NCBI accession to
# 	accession length, taxid, name lineage, and taxid lineage.
# Also maps all lowest-level taxids to organism length, which is the sum of
#  	accession lengths for accessions with that taxid, and lineage info
def get_acc2info(args):
	acc2info, tax2info = {}, {}
	with(open(args.dbinfo, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			rank = get_taxid_rank(taxlin)
			#if rank == 'strain' and acc != 'Unmapped':
			#	taxid += '.1'  # CAMI formatting specification
			#	taxlin += '.1'
			acclen = int(acclen)
			acc2info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in tax2info:
				tax2info[taxid][0] += acclen
			else:
				tax2info[taxid] = [acclen, rank, namelin, taxlin]
	echo('Done reading dbinfo file.', args.verbose)
	return acc2info, tax2info


# Renormalize each taxonomic rank so each rank sums to 100% abundance
def rank_renormalize(tax2abs, only_strains=False):
	rank_totals = {i:0.0 for i in RANKS}  # current rank abundance sums
	mapped_pct = 100.0
	#mapped_pct = 100.0 - (100.0 * tax2abs['Unmapped'][-1])
	for taxid in tax2abs:
		#if taxid == 'Unmapped':
		#	continue
		rank, ab = tax2abs[taxid][1], tax2abs[taxid][-1]
		if only_strains and rank != 'strain':
			continue
		rank_totals[rank] += ab  # add this to the rank sum total

	for taxid in tax2abs:  # here's the normalization
		#if taxid == 'Unmapped':
		#	continue
		rank = tax2abs[taxid][1]
		if only_strains and rank != 'strain':
			continue
		if rank_totals[tax2abs[taxid][1]] != 0.0:
			tax2abs[taxid][-1] /= (rank_totals[tax2abs[taxid][1]] / mapped_pct)
	return tax2abs


# Combines and averages results across all input files,
#  	and packs information into easy-to-write form
def gather_results(args, results):
	echo('Compiling and writing results...', args.verbose)
	rank_results = {i:[] for i in range(len(RANKS))}
	for taxid in results:
		#if 'Unmapped' in taxid:
		#	continue
		ab = results[taxid][-1]  # avg over input files, truncate to 4 digits
		if results[taxid][-1] < 0.00001:
			results[taxid][-1] = '0.00001'
		else:
			results[taxid][-1] = float('%.5f' % results[taxid][-1])
		#results[taxid][-1] = float(math.trunc(
		#	(ab / len(args.infile)) * (10**4))) / (10**4)
		rank = RANKS.index(results[taxid][1])
		if rank == 7:  # strain; add extra CAMI genomeID and OTU fields
			taxid = results[taxid][0]
			cami_genid, cami_otu = taxid, taxid.split('.')[0]
			results[taxid].extend([cami_genid, cami_otu])
		rank_results[rank].append(results[taxid])  # results per rank
	return rank_results


# Writes results out in CAMI format
def write_results(args, rank_results):
	with(open(args.output, 'w')) as outfile:
		# Print some CAMI format header lines
		if args.sampleID == 'NONE':
			outfile.write('@SampleID:' + ','.join(args.infile) + '\n')
		else:
			outfile.write('@SampleID:' + args.sampleID + '\n')
		outfile.write('@Version:MiCoP2-v0.1\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t' +
			'PERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

		for i in range(len(RANKS)):
			lines = rank_results[i]  # all lines to write for this tax level
			# now sort taxids in rank by descending abundance
			if i != 7:  # not strains, which have extra fields
				lines.sort(key=lambda x: 100.0-x[-1])
			else:  # strains
				lines.sort(key=lambda x: 100.0-x[-3])
			if lines == None or len(lines) < 1:
				continue
			for line in lines:
				#if line[4] < 0.0001:  # args.min_abundance:
				#	continue
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')


# Get number of bases hit and total bases using CIGAR string
def parse_cigar(cigar):
	matched_len, total_len, cur = 0, 0, 0
	for ch in cigar:
		if not ch.isalpha():  # ch is a number
			cur = (cur * 10) + int(ch)
		else:  # ch is a letter
			if ch == 'M' or ch == '=':
				matched_len += cur
			total_len += cur
			cur = 0
	return matched_len, total_len


def process_samfile(args):
	contigs_to_accs, contigs_to_bases = {}, {}
	with(open(args.infile, 'r')) as infile:
		for line in infile:
			if line.startswith('@'):
				continue
			splits = line.strip().split('\t')
			cigar = splits[5]
			if cigar == '*':
				continue

			matched_len, total_len = parse_cigar(cigar)
			contig, acc = splits[0], splits[2]
			if acc.startswith('gi|'):
				acc = acc.split('|')[3]
			if '.block' in acc:
				acc = acc.split('.block')[0]
			if contig not in contigs_to_accs:
				contigs_to_accs[contig] = []
				contigs_to_bases[contig] = 0
			contigs_to_accs[contig].append(acc)
			contigs_to_bases[contig] += matched_len
	return contigs_to_accs, contigs_to_bases


# Given all read hits for a read, return their lowest common ancestor (LCA)
def get_lowest_common_ancestor(all_taxids, tax2info):
	lca = ''
	all_taxlins = [tax2info[taxid][-1].split('|') for taxid in all_taxids]
	for i in range(8):
		# get all unique non-blank taxids for this rank
		all_rank_taxids = [taxlin[7 - i] for taxlin in all_taxlins]
		if '' in all_rank_taxids:
			continue
		all_rank_taxids = list(set([taxid for taxid in all_rank_taxids]))
		# if only 1 uniq taxid at this rank, stop and return as LCA
		if len(all_rank_taxids) == 1:
			lca = all_rank_taxids[0]
			break
	return lca


def get_lca_info(lca, tax_to_info):
	if lca in tax_to_info:
		rank, namelin, taxlin = tax_to_info[lca][1:]
		return rank, taxlin, namelin
	rank, taxlin, namelin = '', '', ''
	for tax in tax_to_info:
		tax_namelin, tax_taxlin = tax_to_info[tax][2:]
		# if lca is ancestor of tax, we can get lca information here
		taxlin_splits, namelin_splits = tax_taxlin.split('|'), tax_namelin.split('|')
		if lca not in taxlin_splits:  # lca not a parent of this taxid
			continue
		lca_index = taxlin_splits.index(lca)
		rank = RANKS[lca_index]
		taxlin = '|'.join(taxlin_splits[: lca_index + 1])
		namelin = '|'.join(namelin_splits[: lca_index + 1])
		return rank, taxlin, namelin


def make_lca_to_abs(contigs_to_accs, contigs_to_bases, acc_to_info, tax_to_info):
	lca_to_abs = {}
	for contig in contigs_to_accs:
		acc_list = contigs_to_accs[contig]
		taxids = [acc_to_info[acc][1] for acc in acc_list]
		lca = get_lowest_common_ancestor(taxids, tax_to_info)
		if lca == '':
			continue
		if lca not in lca_to_abs:
			rank, taxlin, namelin = get_lca_info(lca, tax_to_info)
			lca_to_abs[lca] = [lca, rank, taxlin, namelin, 0]
		lca_to_abs[lca][4] += contigs_to_bases[contig]
	return lca_to_abs


def fill_higher_taxa(lca_to_abs):
	tax_to_abs = {}
	for lca in lca_to_abs:
		if lca not in tax_to_abs:
			tax_to_abs[lca] = lca_to_abs[lca]
		else:
			tax_to_abs[lca][4] += lca_to_abs[lca][4]

		# now we also add this abundance to all ancestors of the lca
		taxlin, namelin = lca_to_abs[lca][2:4]
		taxlin_splits, namelin_splits = taxlin.split('|'), namelin.split('|')
		for i in range(len(taxlin_splits)):
			this_taxid = taxlin_splits[i]
			if this_taxid in tax_to_abs:
				tax_to_abs[this_taxid][4] += lca_to_abs[lca][4]
			else:
				this_rank = RANKS[i]
				this_taxlin = '|'.join(taxlin_splits[: i + 1])
				this_namelin = '|'.join(namelin_splits[: i + 1])
				tax_to_abs[this_taxid] = [this_taxid, this_rank,
					this_taxlin, this_namelin, lca_to_abs[lca][4]]
	return tax_to_abs


def main():
	args = parseargs()
	acc_to_info, tax_to_info = get_acc2info(args)
	contigs_to_accs, contigs_to_bases = process_samfile(args)
	lca_to_abs = make_lca_to_abs(contigs_to_accs, contigs_to_bases, acc_to_info, tax_to_info)
	tax_to_abs = fill_higher_taxa(lca_to_abs)
	tax_to_abs = rank_renormalize(tax_to_abs)
	rank_results = gather_results(args, tax_to_abs)
	write_results(args, rank_results)


if __name__ == '__main__':
	main()
#
