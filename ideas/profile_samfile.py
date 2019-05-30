import argparse, math, os, random, subprocess, sys, time
from operator import mul
from functools import reduce


start = time.time()  # start a program timer
__location__ = os.path.realpath(os.path.join(os.getcwd(),
								os.path.dirname(__file__))) + '/'
RANKS = ['superkingdom', 'phylum', 'class', 'order',
		'family', 'genus', 'species', 'strain']


def echo(msg):
	global start
	seconds = time.time() - start
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	hms = "%02d:%02d:%02d" % (h, m, s)
	print('['+hms+'] ' + msg)


def parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Given samfile & profile, extract statistics for TP/FP taxa.')
	parser.add_argument('samfile', help='sam file to process. Required.')
	parser.add_argument('profile', help='CAMI gold standard profile to use for TP/FP. Required.')
	parser.add_argument('--dbinfo', default='AUTO',
		help = 'Location of db_info file. Default: data/subset_db_info.txt')
	parser.add_argument('--no_len_normalization', action='store_true',
		help='Do not normalize abundances by genome length.')
	parser.add_argument('--output', default='samfile_statistics.txt',
		help='Output file name. Default: samfile_statistics.txt')
	parser.add_argument('--pct_id', type=float, default=0.5,
		help='Minimum percent identity from reference to count a hit. Default: 0.5')
	parser.add_argument('--read_cutoff', type=int, default=-1,
		help='Number of reads to count an organism as present.')
	#parser.add_argument('--strain_level', action='store_true',
	#	help='Write output at the strain level as well.')
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
	if args.verbose:
		echo('Reading dbinfo file...')
	acc2info, taxid2info = {}, {}
	with(open(args.dbinfo, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			rank = get_taxid_rank(taxlin)
			if rank == 'strain' and acc != 'Unmapped':
				taxid += '.1'  # CAMI formatting specification
				taxlin += '.1'
			acclen = int(acclen)
			acc2info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in taxid2info:
				taxid2info[taxid][0] += acclen
			else:
				taxid2info[taxid] = [acclen, rank, namelin, taxlin]
	if args.verbose:
		echo('Done reading dbinfo file.')
	return acc2info, taxid2info


def get_pct_id(args, splits):
	cigar = splits[5]  # quality of mapping determined by CIGAR string
	matched_len, total_len, cur = 0, 0, 0
	for ch in cigar:
		if not ch.isalpha():  # ch is a number
			cur = (cur * 10) + int(ch)
		else:  # ch is a letter
			if ch == 'M' or ch == '=':
				matched_len += cur
			total_len += cur
			cur = 0
	return float(matched_len) / float(total_len), total_len


# Parses FLAG field in sam line
def parse_flag(flag, cigar):
	# 1st or 2nd read in a pair, both false if single end
	pair1 = (flag & 1 != 0) and (flag & 64 != 0)
	pair2 = (flag & 1 != 0) and (flag & 128 != 0)
	chimeric = (flag & 2048 != 0)
	# "bad" here means unmapped or cigar string unavailable
	is_bad = (flag & 4 != 0) or (cigar == '*')  # or (flag & 2048 != 0)
	return pair1, pair2, chimeric, is_bad


# for paired end reads, return read hits to references that both paired ends hit
def intersect_read_hits(read_hits, pair1maps, pair2maps):
	if pair1maps == 0 or pair2maps == 0:  # one end unmapped means no intersect
		return [], []
	# gather all ref hits then partition into hits for each pair
	all_ref_hits = [hit[2] for hit in read_hits]
	pair1refs, pair2refs = all_ref_hits[:pair1maps], all_ref_hits[pair1maps:]
	# now intersect the lists and return hits to references in the intersect
	intersect = set([ref for ref in pair1refs if ref in pair2refs])
	#intersect = set([ref for ref in all_ref_hits])
	intersect_hits = [hit for hit in read_hits if hit[2] in intersect]
	return [intersect, intersect_hits]


# Given hits for a read, extract total hit length and quality scores for both
#  	ends of a read pair (if paired), and filter reads using user specification
def clean_read_hits(args, read_hits, pair1maps, pair2maps):
	hitlen, readquals = 0, ''
	filtered_hits = []
	for hit in range(len(read_hits)):
		pct_id, total_len = get_pct_id(args, read_hits[hit])
		read_hits[hit][5] = [pct_id, total_len]  # we no longer need cigar
		# remove user-filtered & chimeric reads, update # of pair end map counts
		if pct_id < args.pct_id or read_hits[hit][1][2]:
			filtered_hits.append(hit)
			if read_hits[hit][1][0]:
				pair1maps -= 1
			elif read_hits[hit][1][1]:
				pair2maps -= 1

		if read_hits[hit][10] != '*':  # first hit for a read / paired end
			readquals += read_hits[hit][10]
			hitlen += len(read_hits[hit][10])
	read_hits = [read_hits[i] for i in range(len(read_hits))
					if i not in filtered_hits]
	return read_hits, [pair1maps, pair2maps, hitlen, readquals]


# Given all hits for a read, decide whether to use it and whether multimapped
# Returns multimapped reads, and if not, taxid and hit length of uniq map
def process_read(args, read_hits, pair1, pair2, pair1maps, pair2maps):
	read_hits, properties = clean_read_hits(args,read_hits,pair1maps,pair2maps)
	pair1maps, pair2maps, hitlen, readquals = properties
	if len(read_hits) == 0:  # all lines filtered
		return [], 'Ambiguous', '', -1
	if pair1 or pair2:  # paired read
		if pair1maps + pair2maps == 1:
			# if uniq mapped to one end and unmapped to other, count mapped end
			return read_hits, read_hits[0][2], readquals, hitlen

		intersect, intersect_hits = intersect_read_hits(
			read_hits, pair1maps, pair2maps)
		if len(intersect) == 0:  # one end unmapped, other multimapped
			return [], 'Ambiguous', '', -1  #we consider this case too ambiguous
		elif len(intersect) == 1:  # read pairs only agree on one taxid
			return intersect_hits, read_hits[0][2], readquals, hitlen  # considered uniq map
		else:  # both ends multimapped, handle later
			return intersect_hits, 'multi', readquals, hitlen

	else:  # single end
		if pair1maps > 1:  # multimapped
			return read_hits, 'multi', readquals, hitlen  # multimapped
		else:
			return False, read_hits[0][2], readquals, hitlen


# Given all read hits for a read, return their lowest common ancestor (LCA)
def get_lowest_common_ancestor(read_hits, taxid2info):
	lca = ''
	all_taxids = [hit[2] for hit in read_hits]
	all_taxlins = [taxid2info[taxid][-1].split('|') for taxid in all_taxids]
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


# Main samfile processing function: parsing and extracting mapping data
def process_samfile(args, acc2info, taxid2info):
	taxids2stats = {}  # taxids to stats and information
	acc2hitpos = {}  # accessions to positions hit, used to calculate coverage
	prev_read, read_hits = '', [] # read tracker, all hits for read (full lines)
	pair1maps, pair2maps = 0, 0  # reads mapped to each pair (single = pair1)
	total_reads = 0

	infile = open(args.samfile, 'r')
	for line in infile:
		if line.startswith('@'):
			continue
		splits = line.strip().split()
		pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
		if chimer or is_bad:  # unmapped or cigar string unavailable
			continue
		splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
		# here we change accession to taxid since we want to
		#  	profile by organisms, not accessions
		acc = splits[2].split('|')[-2]
		splits[2] = acc2info[acc][1]
		splits[7] = acc  # nevermind, save the accession over pnext which we don't use

		read, ref = splits[0], splits[2]
		if read != prev_read:
			total_reads += 1
			if args.verbose and total_reads % 100000 == 0:
				echo('Done processing ' + str(total_reads) + ' read segments.')
			# get uniq hit taxid and hitlen, or multimapped hits intersect
			intersect_hits, taxid, readquals, hitlen = process_read(
				args, read_hits, pair1, pair2, pair1maps, pair2maps)
			# reset these read-specific variables
			prev_read, read_hits, pair1maps, pair2maps = read, [], 0, 0
			if taxid == 'Ambiguous':  # ambiguous mapping (see process_read)
				continue
			lca = get_lowest_common_ancestor(intersect_hits, taxid2info)
			if lca == '':
				lca_index = -1
			else:
				lca_index = taxid2info[intersect_hits[0][2]][-1].split('|').index(lca)
			avg_pctid = sum([hit[5][0] for hit in intersect_hits]) / len(intersect_hits)
			taxids_processed = {}
			for hit in intersect_hits:
				taxlin = taxid2info[hit[2]][-1].split('|')
				namelin = taxid2info[hit[2]][-2].split('|')
				pos, pct_id, hitlen = hit[3], hit[5][0], hit[5][1]

				for i in range(len(taxlin)):
					this_taxid = taxlin[i]
					if this_taxid == '' or this_taxid in taxids_processed:
						continue
					taxids_processed[this_taxid] = ''  # avoid re-processing this taxid
					is_unique = (i <= lca_index)  # unique hit if LCA or an ancestor
					if this_taxid in taxids2stats:
						if is_unique:
							taxids2stats[this_taxid][0] += 1
							taxids2stats[this_taxid][2] += avg_pctid
						else:
							taxids2stats[this_taxid][1] += 1
							taxids2stats[this_taxid][3] += avg_pctid
					else:
						#read to [uniq_hits, multi_hits, uniq_pctid_sum, multi_pctid_sum,
						#  genome length, rank, name lineage, taxid lineage]
						this_rank = RANKS[i]
						this_taxlin = '|'.join(taxlin[:i+1])
						this_namelin = '|'.join(namelin[:i+1])
						if is_unique:
							taxids2stats[this_taxid] = [1, 0, pct_id, 0] + [0,
								this_rank, this_namelin, this_taxlin]
						else:
							taxids2stats[this_taxid] = [0, 1, 0, pct_id] + [0,
								this_rank, this_namelin, this_taxlin]
				# record position and hitlen for lowest level taxid hit
				acc = hit[7]
				if acc not in acc2hitpos:
					acc2hitpos[acc] = [[],[]]
				if taxid != 'multi':
					acc2hitpos[acc][0].append([pos, hitlen])
					#taxids2stats[hit[2]][4].append([pos, hitlen])
				else:
					acc2hitpos[acc][1].append([pos, hitlen])
					#taxids2stats[hit[2]][5].append([pos, hitlen])

		pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
		pair2maps += pair2  # unchanged if pair2 false
		read_hits.append(splits)

	infile.close()
	return taxids2stats, acc2hitpos


# Apply read cutoff and min_abundance thresholds, but latter relative to parent taxa
def prune_tree(args, taxids2stats):
	del_list = []
	for taxid in taxids2stats:
		if taxid == 'Unmapped':
			continue
		taxinfo = taxids2stats[taxid]
		num_reads, taxlin = taxinfo[0], taxinfo[-1]
		num_reads = taxinfo[0]
		if num_reads < args.read_cutoff:
			del_list.append(taxid)
			continue

	# ensure the children of pruned parents are also pruned
	for taxid in taxids2stats:
		taxlin = taxids2stats[taxid][-1].split('|')
		for tax in taxlin:
			if tax in del_list:
				del_list.append(taxid)
	del_list = list(set(del_list))
	for taxid in del_list:
		del taxids2stats[taxid]
	return taxids2stats


# using the read hit positions in acc2hitpos, compute coverage for accesion for
#  	both uniquely and multi mapped reads. then summarize into coverages for taxids.
def compute_coverages(args, taxid2info, taxids2stats, acc2info, acc2hitpos):
	acc2blocks = {}  # collapse hits into blocks of coverage
	for acc in acc2hitpos:
		acc2blocks[acc] = [[],[]]
		for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
			hits = acc2hitpos[acc][j]
			if len(hits) == 0:
				continue
			hits.sort(key = lambda x: x[0])  # sort by increasing position
			for hit_start, hit_len in hits:
				hit_start, hit_len = int(hit_start), int(hit_len)
				hit_end = hit_start + hit_len
				# three cases: hit creates new block, extends it, or doesn't extend it
				if len(acc2blocks[acc][j]) == 0:
					acc2blocks[acc][j].append([hit_start, hit_end])
				else:
					blocks_end = acc2blocks[acc][j][-1][1]  # end of last covered block
					if hit_start > blocks_end:
						acc2blocks[acc][j].append([hit_start, hit_end])
					elif hit_end > blocks_end:  # hit starts in block but extends it
						acc2blocks[acc][j][-1][1] = hit_end
					# last case is that this hit is subsumed by last block -- do nothing

	acc2coverage = {}
	for acc in acc2blocks:
		acc2coverage[acc] = [[],[]]
		for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
			if len(acc2blocks[acc][j]) == 0:
				continue
			block_lengths = [block[1] - block[0] for block in acc2blocks[acc][j]]
			covered_bases = sum(block_lengths)
			total_bases = acc2info[acc][0]
			coverage = float(covered_bases) / float(total_bases)
			acc2coverage[acc][j] = [coverage, covered_bases, total_bases, block_lengths]

	taxids2coverage = {}
	for taxid in taxid2info:
		if taxid not in taxids2stats:  # a taxid that wasn't mapped to
			continue
		# accs hit, total accs, bases covered, total bases in accs hit & in all accs, block lengths
		taxids2coverage[taxid] = [[0, 0, 0, 0, 0, []], [0, 0, 0, 0, 0, []]]

		taxid_len = taxid2info[taxid][0]
		taxids2coverage[taxid][0][4] = taxid_len
		taxids2coverage[taxid][1][4] = taxid_len
		for acc in acc2coverage:
			if acc2info[acc][1] != taxid:
				continue
			for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
				taxids2coverage[taxid][j][1] += 1
				#taxids2coverage[taxid][j][4] += acc2info[acc][0]
				if acc in acc2coverage and len(acc2coverage[acc][j]) != 0:
					taxids2coverage[taxid][j][0] += 1
					taxids2coverage[taxid][j][2] += acc2coverage[acc][j][1]
					taxids2coverage[taxid][j][3] += acc2coverage[acc][j][2]
					taxids2coverage[taxid][j][5].extend(acc2coverage[acc][j][3])
	return acc2coverage, taxids2coverage


# given a gold stanard profile, compute whether taxids are true or false positives
def compute_tpfp(args, taxids2stats, acc2covg, acc2info):
	tp_taxids = {}  # true positive taxIDs extracted from profile
	with(open(args.profile, 'r')) as infile:
		for line in infile:
			if line.startswith('@') or line.startswith('#') or len(line) < 5:
				continue
			taxlin = line.split('\t')[2].split('|')
			for taxid in taxlin:
				tp_taxids[taxid] = True

	acc2tpfp, taxids2tpfp = {}, {}
	for taxid in taxids2stats:
		if taxid in tp_taxids:
			taxids2tpfp[taxid] = True
		else:
			taxids2tpfp[taxid] = False
	for acc in acc2covg:
		acc_taxid = acc2info[acc][1]
		if acc_taxid in tp_taxids:
			acc2tpfp[acc] = True
		else:
			acc2tpfp[acc] = False
	return acc2tpfp, taxids2tpfp


# compute given quantile values for the stats for taxids in taxids2stats,
#  	using either only true positives or false positives
def compute_stats_quantiles(quantiles, taxids2stats, taxids2tpfp, true_positives):
	if true_positives:
		taxstats = {k: v for k,v in taxids2stats.items() if taxids2tpfp[k]}
	else:
		taxstats = {k: v for k,v in taxids2stats.items() if not taxids2tpfp[k]}
	# taxid to [uniq_hits, multi_hits, uniq_pctid_sum, multi_pctid_sum, taxonomy info]

	# gather and sort this info for tp or fp taxa
	all_uniq_hits, all_multi_hits, all_uniq_pctid, all_multi_pctid = [], [], [], []
	for taxid in taxstats:
		uniq_hits, multi_hits, uniq_pctid_sum, multi_pctid_sum = taxstats[taxid][:4]
		if uniq_hits > 0:
			uniq_avg_pctid = uniq_pctid_sum / float(uniq_hits)
		else:
			uniq_avg_pctid = 0.0
		if multi_hits > 0:
			multi_avg_pctid = multi_pctid_sum / float(multi_hits)
		else:
			multi_avg_pctid = 0.0
		all_uniq_hits.append(uniq_hits)
		all_multi_hits.append(multi_hits)
		all_uniq_pctid.append(uniq_avg_pctid)
		all_multi_pctid.append(multi_avg_pctid)
	all_uniq_pctid = [float('%.5f'%(k)) for k in all_uniq_pctid]
	all_multi_pctid = [float('%.5f'%(k)) for k in all_multi_pctid]
	all_uniq_hits.sort() ; all_multi_hits.sort() ; all_uniq_pctid.sort() ; all_multi_pctid.sort()

	# compute quantile values: get indices for each quantile, extract info into dict
	stats_quantiles = {}
	quant_ind = [int(quantiles[i] * len(all_uniq_hits)) for i in range(len(quantiles))]
	stats_quantiles['Unique hits'] = [all_uniq_hits[ind] for ind in quant_ind]
	stats_quantiles['Multi hits'] = [all_multi_hits[ind] for ind in quant_ind]
	stats_quantiles['Unique avg. pctid'] = [all_uniq_pctid[ind] for ind in quant_ind]
	stats_quantiles['Multi avg. pctid'] = [all_multi_pctid[ind] for ind in quant_ind]
	return stats_quantiles


# Compute given quantile values for coverage statistics of TP/FP TaxIDs
def compute_taxid_covg_quantiles(quantiles, taxids2covg, taxids2tpfp, true_positives):
	if true_positives:
		tax_covg = {k: v for k,v in taxids2covg.items() if taxids2tpfp[k]}
	else:
		tax_covg = {k: v for k,v in taxids2covg.items() if not taxids2tpfp[k]}
	# taxid to [accs hit, total accs, bases covered, bases in accs hit & total, block_lengths]

	if len(tax_covg) == 0:
		sys.exit('Error: one or more classes is empty.')
	all_covg_in_hits, all_covg_overall, all_block_hits, all_block_lens = [[],[]],[[],[]],[[],[]],[[],[]]
	for taxid in tax_covg:
		for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
			if len(tax_covg[taxid][j]) == 0 or tax_covg[taxid][j][0] == 0:
				continue  # no hit accs for this mapping status for this taxid
			hit_accs, total_accs, cov_bases, hit_acc_bases, total_bases, block_lens = tax_covg[taxid][j]
			pct_accs_hit = float(hit_accs) / float(total_accs)
			# compute coverage % in only the accesions hit by reads, and overall bases covered
			covg_in_hits = float(cov_bases) / float(hit_acc_bases)
			covg_overall = float(cov_bases) / float(total_bases)
			block_hits = len(block_lens)
			all_covg_in_hits[j].append(covg_in_hits)
			all_covg_overall[j].append(covg_overall)
			all_block_hits[j].append(block_hits)
			all_block_lens[j].extend(block_lens)
	all_covg_in_hits[0] = [float('%.5f'%(k)) for k in all_covg_in_hits[0]]
	all_covg_in_hits[1] = [float('%.5f'%(k)) for k in all_covg_in_hits[1]]
	all_covg_overall[0] = [float('%.5f'%(k)) for k in all_covg_in_hits[0]]
	all_covg_overall[1] = [float('%.5f'%(k)) for k in all_covg_in_hits[1]]
	all_covg_in_hits[0].sort() ; all_covg_overall[0].sort()
	all_block_hits[0].sort() ; all_block_lens[0].sort()
	all_covg_in_hits[1].sort() ; all_covg_overall[1].sort()
	all_block_hits[1].sort() ; all_block_lens[1].sort()

	# compute quantile values: get indices for each quantile, extract info into dict
	covg_quantiles = {}
	quant_ind_uniq = [int(quantiles[i] * len(all_covg_in_hits[0])) for i in range(len(quantiles))]
	quant_ind_multi = [int(quantiles[i] * len(all_covg_in_hits[1])) for i in range(len(quantiles))]
	block_quant_ind_uniq = [int(quantiles[i] * len(all_block_lens[0])) for i in range(len(quantiles))]
	block_quant_ind_multi = [int(quantiles[i] * len(all_block_lens[1])) for i in range(len(quantiles))]
	covg_quantiles['TaxID uniquely mapped reads coverage over hit accessions'] = [
		all_covg_in_hits[0][ind] for ind in quant_ind_uniq]
	covg_quantiles['TaxID uniquely mapped reads coverage over all accessions'] = [
		all_covg_overall[0][ind] for ind in quant_ind_uniq]
	covg_quantiles['TaxID uniquely mapped reads number of blocks'] = [
		all_block_hits[0][ind] for ind in quant_ind_uniq]
	covg_quantiles['TaxID uniquely mapped reads block lengths'] = [
		all_block_lens[0][ind] for ind in block_quant_ind_uniq]
	covg_quantiles['TaxID multimapped reads coverage over hit accessions'] = [
		all_covg_in_hits[1][ind] for ind in quant_ind_multi]
	covg_quantiles['TaxID multimapped reads coverage over all accessions'] = [
		all_covg_overall[1][ind] for ind in quant_ind_multi]
	covg_quantiles['TaxID multimapped reads number of blocks'] = [
		all_block_hits[1][ind] for ind in quant_ind_multi]
	covg_quantiles['TaxID multimapped reads block lengths'] = [
		all_block_lens[1][ind] for ind in block_quant_ind_multi]
	return covg_quantiles


# Compute given quantile values for coverage statistics of TP/FP accessions
def compute_acc_covg_quantiles(quantiles, acc2covg, acc2tpfp, true_positives):
	if true_positives:
		acc_covg = {k: v for k,v in acc2covg.items() if acc2tpfp[k]}
	else:
		acc_covg = {k: v for k,v in acc2covg.items() if not acc2tpfp[k]}
	# acc_covg has [coverage_pct, covered_bases, total_bases, block_lengths]

	if len(acc_covg) == 0:
		sys.exit('Error: one or more classes is empty.')
	all_coverage, all_block_hits, all_block_lens = [[],[]], [[],[]], [[],[]]
	for acc in acc_covg:
		for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
			if len(acc_covg[acc][j]) == 0 or acc_covg[acc][j][0] == 0.0:
				continue  # no coverage for this accession for this mapping status
			coverage, cov_bases, total_bases, block_lengths = acc_covg[acc][j]
			block_hits = len(block_lengths)
			all_coverage[j].append(coverage)
			all_block_hits[j].append(block_hits)
			all_block_lens[j].extend(block_lengths)
	all_coverage[0] = [float('%.5f'%(k)) for k in all_coverage[0]]
	all_coverage[1] = [float('%.5f'%(k)) for k in all_coverage[1]]
	all_coverage[0].sort() ; all_block_hits[0].sort() ; all_block_lens[0].sort()
	all_coverage[1].sort() ; all_block_hits[1].sort() ; all_block_lens[1].sort()

	# compute quantile values: get indices for each quantile, extract info into dict
	covg_quantiles = {}
	quant_ind_uniq = [int(quantiles[i] * len(all_coverage[0])) for i in range(len(quantiles))]
	quant_ind_multi = [int(quantiles[i] * len(all_coverage[1])) for i in range(len(quantiles))]
	block_quant_ind_uniq = [int(quantiles[i] * len(all_block_lens[0])) for i in range(len(quantiles))]
	block_quant_ind_multi = [int(quantiles[i] * len(all_block_lens[1])) for i in range(len(quantiles))]
	covg_quantiles['Accesion uniquely mapped reads coverage'] = [
		all_coverage[0][ind] for ind in quant_ind_uniq]
	covg_quantiles['Accession uniquely mapped reads number of blocks'] = [
		all_block_hits[0][ind] for ind in quant_ind_uniq]
	covg_quantiles['Accession uniquely mapped reads block lengths'] = [
		all_block_lens[0][ind] for ind in block_quant_ind_uniq]
	covg_quantiles['Accesion multimapped reads coverage'] = [
		all_coverage[1][ind] for ind in quant_ind_multi]
	covg_quantiles['Accession multimapped reads number of blocks'] = [
		all_block_hits[1][ind] for ind in quant_ind_multi]
	covg_quantiles['Accession multimapped reads block lengths'] = [
		all_block_lens[1][ind] for ind in block_quant_ind_multi]
	return covg_quantiles


# write the information out in a human-readable / user-friendly format
def write_results(args, taxids2stats, taxids2covg, taxids2tpfp, acc2covg, acc2tpfp, acc2info):
	# for given quantiles, compute quantiles for all these metrics for both TPs and FPs
	prevent_oob = 0.999999999  # prevent index out of bounds error for list.index(len(list)*1)
	#quantiles = [0.0, 0.25, 0.5, 0.75, prevent_oob]
	quantiles = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, prevent_oob]
	true_stats_quantiles = compute_stats_quantiles(quantiles, taxids2stats, taxids2tpfp, True)
	false_stats_quantiles = compute_stats_quantiles(quantiles, taxids2stats, taxids2tpfp, False)
	true_taxid_covg_quantiles = compute_taxid_covg_quantiles(quantiles, taxids2covg, taxids2tpfp, True)
	false_taxid_covg_quantiles = compute_taxid_covg_quantiles(quantiles, taxids2covg, taxids2tpfp, False)
	true_acc_covg_quantiles = compute_acc_covg_quantiles(quantiles, acc2covg, acc2tpfp, True)
	false_acc_covg_quantiles = compute_acc_covg_quantiles(quantiles, acc2covg, acc2tpfp, False)

	with(open(args.output, 'w')) as outfile:
		if quantiles[0] == 0.0: quantiles[0] = 'min'  # for printing purposes
		if quantiles[-1] == prevent_oob: quantiles[-1] = 'max'  # for printing purposes
		# write the quantile results computed above
		outfile.write('True positive results for quantiles: ' + str(quantiles) + '\n')
		for metric in true_stats_quantiles:
			outfile.write(metric + ': ' + str(true_stats_quantiles[metric]) + '\n')
		for metric in true_taxid_covg_quantiles:
			outfile.write(metric + ': ' + str(true_taxid_covg_quantiles[metric]) + '\n')
		for metric in true_acc_covg_quantiles:
			outfile.write(metric + ': ' + str(true_acc_covg_quantiles[metric]) + '\n')
		outfile.write('\nFalse positive results for quantiles: ' + str(quantiles) + '\n')
		for metric in false_stats_quantiles:
			outfile.write(metric + ': ' + str(false_stats_quantiles[metric]) + '\n')
		for metric in false_taxid_covg_quantiles:
			outfile.write(metric + ': ' + str(false_taxid_covg_quantiles[metric]) + '\n')
		for metric in false_acc_covg_quantiles:
			outfile.write(metric + ': ' + str(false_acc_covg_quantiles[metric]) + '\n')

		# write the per-taxid information
		outfile.write('\n\n\nPer taxid information: \n\n')
		rank_ordered_stats = [[] for i in range(len(RANKS))]
		for i in range(len(RANKS)):
			#if i == 7 and not args.strain_level:
			#	continue
			for taxid in taxids2stats:
				rank, namelin, taxlin = taxids2stats[taxid][5:]
				if rank != RANKS[i]:
					continue
				taxstats = taxids2stats[taxid][:4]
				# convert pct id sums to averages
				if taxstats[0] != 0:
					taxstats[2] = taxstats[2] / float(taxstats[0])
				if taxstats[1] != 0:
					taxstats[3] = taxstats[3] / float(taxstats[1])
				if taxids2tpfp[taxid]:
					status = 'True Positive'
				else:
					status = 'False Positive'
				# [taxid, rank, status, taxlin, namelin, uniq hits, multi hits, uniq pctid, multi pctid]
				rank_ordered_stats[i].append([taxid, rank, status, taxlin, namelin] + taxstats)

		for i in range(len(rank_ordered_stats)):
			rank_ordered_stats[i].sort(key = lambda x: x[5], reverse=True)  # sort by uniq reads hit
			for line in rank_ordered_stats[i]:
				line = [str(k) for k in line]
				outfile.write('\t'.join(line[:5]) + '\n')
				outfile.write('[uniq hits, multi hits, uniq pctid, multi pctid] == ')
				outfile.write('\t'.join(line[5:]) + '\n')
				if line[0] in taxids2covg:
					for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
						if j == 0:
							outfile.write('Uniquely mapped reads ')
						else:
							outfile.write('Multimapped reads ')
						if len(taxids2covg[line[0]][j]) == 0 or taxids2covg[line[0]][j][0] == 0:
							pct_accs_hit, covg_in_hits, covg_overall, block_hits = 0.0, 0.0, 0.0, 0
						else:
							x = taxids2covg[line[0]][j]
							hit_accs, tot_accs, cov_bases, hit_acc_bases, tot_acc_bases, block_lens = x
							pct_accs_hit = float(hit_accs) / float(tot_accs)
							# compute coverage % in only the accs hit by reads, and overall bases covered
							covg_in_hits = float(cov_bases) / float(hit_acc_bases)
							covg_overall = float(cov_bases) / float(tot_acc_bases)
							block_hits = len(block_lens)
							avg_block_hit = float(sum(block_lens)) / float(len(block_lens))
						to_write = str([str(k) for k in
							[pct_accs_hit, covg_in_hits, covg_overall, block_hits, avg_block_hit]])
						outfile.write('[Pct. accs hit, Avg. covg. in accs hit, ')
						outfile.write('Avg. covg. over all accs, Num. blocks hit, ')
						outfile.write('Avg. block hit len.] == \n' + to_write + '\n')
				outfile.write('\n\n')

		# write the per-accession information
		outfile.write('\n\n\nPer accession information: \n\n')
		acc_covgs = [[acc] + acc2covg[acc] for acc in acc2covg]
		acc_covgs.sort(key = lambda x: x[1], reverse=True)  # sort by coverage percentage
		for acc_line in acc_covgs:
			accession = acc_line[0]
			if acc2tpfp[accession]:
				status = 'True positive'
			else:
				status = 'False positive'
			acclen, taxid, namelin, taxlin = acc2info[accession]
			outfile.write('\t'.join([accession, taxid, status, taxlin, namelin]) + '\n')
			# acc_covgs has [coverage_pct, covered_bases, total_bases, block_lengths]
			for j in range(2):  # j=0 --> uniq mapped, j=1 --> multi mapped
				if j == 0:
					outfile.write('Uniquely mapped reads ')
				else:
					outfile.write('Multimapped reads ')
				if len(acc2covg[accession][j]) == 0:
					covg_pct, cov_bases, total_bases, block_lengths = 0.0, 0, 0, [0]
				else:
					covg_pct, cov_bases, total_bases, block_lengths = acc2covg[accession][j]
				num_blocks = len(block_lengths)
				avg_block_len = float(sum(block_lengths)) / float(len(block_lengths))
				outfile.write('[Coverage pct., Num. blocks, Avg. block len.] == \n')
				outfile.write(str([str(k) for k in [covg_pct, num_blocks, avg_block_len]]) + '\n')
			outfile.write('\n')


def main(args = None):
	if args == None:
		args = profile_parseargs()
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		print('Error: --pct_id must be between 0.0 and 1.0, inclusive.')
		sys.exit()
	if args.dbinfo == 'AUTO':
		args.dbinfo = __location__ + 'data/db_info.txt'
	open(args.output, 'w').close()  # test to see if writeable

	# maps NCBI accession to length, taxid, name lineage, taxid lineage
	acc2info, taxid2info = get_acc2info(args)
	taxids2stats, acc2hitpos = process_samfile(args, acc2info, taxid2info)
	taxids2stats = prune_tree(args, taxids2stats)

	acc2covg, taxids2covg = compute_coverages(args, taxid2info, taxids2stats, acc2info, acc2hitpos)
	acc2tpfp, taxids2tpfp = compute_tpfp(args, taxids2stats, acc2covg, acc2info)
	write_results(args, taxids2stats, taxids2covg, taxids2tpfp, acc2covg, acc2tpfp, acc2info)


if __name__ == '__main__':
	args = parseargs()
	main(args)
#
