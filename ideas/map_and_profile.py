import argparse, ast, math, os, subprocess, sys, time
import numpy as np
from scipy.stats import kstest


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


def profile_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Compute abundance estimations for species in a sample.')
	parser.add_argument('infiles', nargs='+',
		help='sam or reads file(s) (space-delimited if multiple). Required.')
	parser.add_argument('--db', default='NONE',
		help='Path to database from select_db.py. Required if read files given.')
	parser.add_argument('--dbinfo', default='AUTO',
		help = 'Location of db_info file. Default: data/subset_db_info.txt')
	#parser.add_argument('--min_abundance', type=float, default=10**-4,
	#	help='Minimum abundance for a taxa to be included in the results.')
	parser.add_argument('--length_normalization', action='store_true',
		help='Normalize abundances by genome length.')
	parser.add_argument('--no_rank_renormalization', action='store_true',
		help='Do not renormalize abundances to 100 percent at each rank,\
				for instance if an organism has a species but not genus label.')
	parser.add_argument('--output', default='abundances.tsv',
		help='Output abundances file. Default: abundances.tsv')
	parser.add_argument('--pct_id', type=float, default=0.5,
		help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--no_quantify_unmapped', action='store_false',
		help='Factor in unmapped reads in abundance estimation.')
	parser.add_argument('--read_cutoff', type=int, default=100,
		help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--no_strain_level', action='store_true',
		help='Write output at the strain level as well.')
	parser.add_argument('--uniq_blocks_file',
		help = 'File mapping accesions to unique blocks. Required.')
	#parser.add_argument('--uniq_covg_cutoff', default = 0.2, type = float,
	#	help='Coverage pct. by uniquely-mapped reads needed to count as present.')
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


# Combine the unique hit positions into blocks of contiguous unique hits
def read_unique_blocks(uniq_blocks_file, acc2info):
	acc2uniqblocks = {}
	with(open(uniq_blocks_file, 'r')) as infile:
		for line in infile:
			splits = line.strip().split('\t')
			if len(splits) == 1:
				acc2uniqblocks[splits[0]] = []
			else:
				blocks = [ast.literal_eval(i) for i in splits[1:]]
				acc2uniqblocks[splits[0]] = blocks
	return acc2uniqblocks


# Get pct match of read hit using CIGAR string
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
	return float(matched_len) / float(total_len), matched_len


# Parses FLAG field in sam line
def parse_flag(flag, cigar):
	# 1st or 2nd read in a pair, both false if single end
	pair1 = (flag & 1 != 0) and (flag & 64 != 0)
	pair2 = (flag & 1 != 0) and (flag & 128 != 0)
	chimeric = (flag & 2048 != 0)
	# "bad" here means unmapped or cigar string unavailable
	is_bad = (flag & 4 != 0) or (cigar == '*')  # or (flag & 2048 != 0)
	return pair1, pair2, chimeric, is_bad


# Given hits for a read, extract total hit length and quality scores for both
#  	ends of a read pair (if paired), and filter reads using user specification
def filter_read_hits(args, read_hits, pair1maps, pair2maps):
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

	read_hits = [read_hits[i] for i in range(len(read_hits))
					if i not in filtered_hits]
	return read_hits, pair1maps, pair2maps


# Given all hits for a read, decide whether to use it and whether multimapped
# Returns multimapped reads, and if not, taxid and hit length of uniq map
def intersect_read_hits(args, read_hits, pair1, pair2, pair1maps, pair2maps):
	if len(read_hits) == 0:  # all lines filtered
		return [], 'Ambiguous'
	if pair1 or pair2:  # paired read
		if pair1maps + pair2maps == 1:
			# if uniq mapped to one end and unmapped to other, count mapped end
			return read_hits, read_hits[0][7]

		# for paired end reads, find read hits to references that both paired ends hit
		if pair1maps == 0 or pair2maps == 0:  # one end unmapped means no intersect
			return [], 'Ambiguous'
		# gather all ref hits then partition into hits for each pair
		all_ref_hits = [hit[7] for hit in read_hits]
		pair1refs, pair2refs = all_ref_hits[:pair1maps], all_ref_hits[pair1maps:]
		# now intersect the lists and return hits to references in the intersect
		intersect = set([ref for ref in pair1refs if ref in pair2refs])
		#intersect = set([ref for ref in all_ref_hits])
		intersect_hits = [hit for hit in read_hits if hit[7] in intersect]

		if len(intersect) == 0:  # one end unmapped, other multimapped
			return [], 'Ambiguous'  #we consider this case too ambiguous
		elif len(intersect) == 1:  # read pairs only agree on one taxid
			return intersect_hits, read_hits[0][7]  # considered uniq map
		else:  # both ends multimapped, handle later
			return intersect_hits, 'multi'

	else:  # single end
		if pair1maps > 1:  # multimapped
			return read_hits, 'multi'  # multimapped
		else:
			return False, read_hits[0][7]


# Given all read hits for a read, return their lowest common ancestor (LCA)
def get_lowest_common_ancestor(read_hits, tax2info):
	lca = ''
	all_taxids = [hit[7] for hit in read_hits]
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


# Given hit information and LCA, update abundance info for relevant taxa
def process_hit(hit, taxid, lca_index, avg_hitlen, acc2hitpos, tax2abs, tax2info):
	taxids_processed = {}
	# record position and hitlen for lowest level taxid hit
	acc, pos, pct_id, hitlen = hit[2], int(hit[3]), hit[5][0], hit[5][1]
	if acc not in acc2hitpos:
		acc2hitpos[acc] = [[],[]]
	if taxid != 'multi':
		acc2hitpos[acc][0].append([pos, hitlen])
	else:
		acc2hitpos[acc][1].append([pos, hitlen])

	taxlin = tax2info[hit[7]][-1].split('|')
	namelin = tax2info[hit[7]][-2].split('|')
	for i in range(len(taxlin)):
		this_taxid = taxlin[i]
		if this_taxid == '' or this_taxid in taxids_processed:
			continue
		taxids_processed[this_taxid] = ''  # avoid re-processing this taxid
		is_unique = (i <= lca_index)  # unique hit if LCA or an ancestor
		if this_taxid in tax2abs:
			if this_taxid in tax2info:  # lowest level taxid
				if is_unique:
					tax2abs[this_taxid][0] += 1
					tax2abs[this_taxid][1] += avg_hitlen
			else:
				if is_unique:
					tax2abs[this_taxid][0] += 1
					tax2abs[this_taxid][1] += avg_hitlen
					if taxid == 'multi':
						if acc not in tax2abs[this_taxid][2]:
							tax2abs[this_taxid][2][acc] = [len(acc2hitpos[acc][1]) - 1]
						else:
							tax2abs[this_taxid][2][acc].append(len(acc2hitpos[acc][1]) - 1)
			#else: tax2abs[this_taxid][1] += 1
		else:
			#read to [uniq_hits, uniq_hit_bases,
			#  3rd field (see below), rank, name lineage, taxid lineage]
			if this_taxid in tax2info:  # lowest level taxid
				# 3rd field is accession length -- for covg purposes
				if is_unique:
					tax2abs[this_taxid] = [1, avg_hitlen] + tax2info[this_taxid]
				#else:
				#	tax2abs[this_taxid] = [0, 0] + tax2info[this_taxid]
			else:
				# 3rd field is dict with counts of reads unique to this
				#  taxa but multimapped to children, used as proxy
				#  evidence of unique read "coverage" for higher taxa
				this_rank = RANKS[i]
				this_taxlin = '|'.join(taxlin[:i+1])
				this_namelin = '|'.join(namelin[:i+1])
				if is_unique:
					tax2abs[this_taxid] = [1, avg_hitlen, {},
						this_rank, this_namelin, this_taxlin]
					if taxid == 'multi':
						tax2abs[this_taxid][2][acc] = [len(acc2hitpos[acc][1]) - 1]
				#else:
				#	tax2abs[this_taxid] = [0, 0, {},
				#		this_rank, this_namelin, this_taxlin]
	return tax2abs, taxids_processed, acc2hitpos


# Given the possible read hits & info for this read, update accumulated global results
def update_global_hits(args, hits, taxid, global_results, tax2info, acc2info):
	tax2abs, multimapped, acc2hitpos = global_results
	if taxid == 'Ambiguous':  # ambiguous mapping (see intersect_read_hits)
		if not args.no_quantify_unmapped:
			if 'Unmapped' in tax2abs:
				tax2abs['Unmapped'][0] += 1.0  # reads hit
			else:  # also store lineage for taxid
				tax2abs['Unmapped'] = ([1.0, 0.0] +
					tax2info['Unmapped'])
		return tax2abs, multimapped, acc2hitpos
	if taxid != 'multi' and args.length_normalization:
		for i in range(len(hits)):
			hits[i][5][1] /= tax2info[hits[i][7]][0] # normalize by genome length
	lca = get_lowest_common_ancestor(hits, tax2info)
	if lca == '':
		lca_index = -1
	else:
		lca_index = tax2info[hits[0][7]][-1].split('|').index(lca)
	avg_hitlen = float(sum([hit[5][1] for hit in hits])) / len(hits)
	for hit in hits:
		tax2abs, taxids_processed, acc2hitpos = process_hit(hit, taxid, lca_index,
			avg_hitlen, acc2hitpos, tax2abs, tax2info)

	if taxid == 'multi':  # store taxids hit and length (in bases) of hits
		#hits = [[hit[7], len(hit[9])] for hit in hits]
		hits = [avg_hitlen] + [hit[7] for hit in hits]
		multimapped.append(hits)
	return tax2abs, multimapped, acc2hitpos


# Runs minimap2, processes output, and returns abundance information for taxids,
#  and multimapping and coverage information for post-processing.
# This is the general iterator function; it calls subroutines to process
#  individual lines of input
def map_and_process(args, infile, acc2info, tax2info):
	# fields tracked through lines: current read, previous read, hits to current
	#    read, hit counts for first and second paired ends, total reads
	read, prev_read, read_hits, pair1maps, pair2maps, total_reads = '', '', [], 0, 0, 0
	# dict of taxids to abundances, multimapped reads, read hit positions for accesions
	tax2abs, multimapped, acc2hitpos = {}, [], {}

	samfile = False  # whether reading from existing sam file
	if infile.endswith('sam'):  # input stream from sam file
		samfile = True
		instream = open(infile, 'r')
	else:  # run minimap2 and stream its output as input
		mapper = subprocess.Popen([__location__ + 'minimap2/minimap2', '-ax',
			'sr', '-t', '4', '-2', '-n' '1', '--secondary=yes',
			args.db, infile], stdout=subprocess.PIPE, bufsize=1)
		instream = iter(mapper.stdout.readline, "")

	for line in instream:
		if not samfile:
			line = line.decode('utf-8')
			if not line:  # process finished
				break

		if line.startswith('@'):
			if not args.no_quantify_unmapped:
				if 'Unmapped' in tax2abs:
					tax2abs['Unmapped'][0] += 1.0  # reads hit
				else:  # also store lineage for taxid
					tax2abs['Unmapped'] = [1.0, 0.0] + tax2info['Unmapped']
			continue
		splits = line.strip().split()
		pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
		if is_bad:  # unmapped or cigar string unavailable
			continue
		splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
		if 'gi|' in splits[2]:
			acc = splits[2].split('|')[3]
			splits[2] = acc
		else:
			acc = splits[2]
		splits[7] = acc2info[acc][1]  # replace PNEXT field (unused) w/ taxid

		read = splits[0]
		if read != prev_read:
			total_reads += 1
			if total_reads % 100000 == 0:
				echo('Done processing ' + str(total_reads) + ' reads.', args.verbose)
			# get uniq hit taxid or multimapped hits intersect
			read_hits, pair1maps, pair2maps = filter_read_hits(
				args, read_hits, pair1maps, pair2maps)
			intersect_hits, taxid = intersect_read_hits(
				args, read_hits, pair1, pair2, pair1maps, pair2maps)
			# update accumulated global results using the read hits information
			global_results = [tax2abs, multimapped, acc2hitpos]
			tax2abs, multimapped, acc2hitpos = update_global_hits(
				args, intersect_hits, taxid, global_results, tax2info, acc2info)
			# reset these read-specific variables
			prev_read, read_hits, pair1maps, pair2maps = read, [], 0, 0

		pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
		pair2maps += pair2  # unchanged if pair2 false
		read_hits.append(splits)

	if samfile:
		instream.close()
	else:
		mapper.stdout.close()
		mapper.wait()
	# store percentage of unmapped reads
	if not args.no_quantify_unmapped:
		tax2abs['Unmapped'][1] = tax2abs['Unmapped'][0] / float(total_reads)
	return tax2abs, multimapped, acc2hitpos


# Given the unique blocks for each accession, break the accessions into their
#    unique blocks, create new names and acc2info entries for those blocks, and
#    assign the appropriate hits to them, with new (relative) positions
def break_accs_into_uniq_blocks(acc2hitpos, acc2info, tax2info, tax2abs, acc2uniqblocks):
	# accs to delete and add to acc2hitpos, and new taxid lengths using just uniq blocks
	del_list, add_list, taxid_lens = {}, {}, {}
	for acc in acc2hitpos:
		taxid = acc2info[acc][1]
		if taxid not in taxid_lens:
			taxid_lens[taxid] = 0
		if acc not in acc2uniqblocks:
			del_list[acc] = True
			continue
		block_num = 1  # unique block number for this accession
		hitpos, hitlen = -1, -1
		uniq_blks = acc2uniqblocks[acc]
		for bl_start, bl_end in uniq_blks:
			block_acc_name = acc + '.block' + str(block_num)
			hit_list = []
			while len(acc2hitpos[acc][0]) > 0:
				hitend = hitpos + hitlen
				if hitpos > bl_end:
					break
				else:
					if hitend >= bl_start:
						hitpos -= bl_start  # get position relative to block
						hit_list.append([hitpos, hitlen])
					hitpos, hitlen = acc2hitpos[acc][0].pop(0)
			#acc2hitpos[block_acc_name] = hit_list
			add_list[block_acc_name] = hit_list
			acc2info[block_acc_name] = acc2info[acc][0:]
			bl_size = bl_end - bl_start
			acc2info[block_acc_name][0] = bl_size  # set block size
			taxid_lens[taxid] += bl_size
			block_num += 1
		del_list[acc] = True

	# update acc and taxid information using just unique blocks
	for acc in del_list:
		del acc2hitpos[acc]  # remove old entry so it isn't used
		del acc2info[acc]
	for acc in add_list:
		acc2hitpos[acc] = [add_list[acc], [[0,1]]]  # uniq hits and multi hits (none in this case)
	for taxid in tax2abs:
		if taxid not in tax2info:  # not leaf node
			continue
		#if taxid not in taxid_lens:
		#	tax2info[taxid][0] = 0
		#	tax2abs[taxid][4] = 0
		#else:
		tax2info[taxid][0] = taxid_lens[taxid]
		tax2abs[taxid][2] = taxid_lens[taxid]
	return acc2hitpos, acc2info, tax2info, tax2abs


# Divide abundances of multimapped reads according to uniquely mapped reads
#  	portion for each of the hit organisms
def resolve_multi_prop(args, tax2abs, multimapped, tax2info):
	# remove hits to pruned nodes, and exit if no multimapped reads remain
	for i in range(len(multimapped)):
		multimapped[i] = [multimapped[i][0]] + [
			hit for hit in multimapped[i][1:] if hit in tax2abs]
	multimapped = [read for read in multimapped if len(read) > 1]
	if len(multimapped) == 0:
		return tax2abs

	to_add = {}  # abundances to add to total after all reads processed
	for read_hits in multimapped:
		# get abundances of all taxids in read_hits
		#all_taxids = list(set([hit[0] for hit in read_hits
		#	if hit[0] in tax2abs]))
		all_taxids = list(set(read_hits[1:]))
		if len(all_taxids) == 0:  # all hits were to taxids with no unique hits
			continue
		taxid_abs = [tax2abs[tax][1] for tax in all_taxids]

		## now get cumulative proportional abs. of taxids relative to each other
		sumabs = sum(taxid_abs)
		if sumabs == 0.0:
			continue
		proportions = [ab / sumabs for ab in taxid_abs]
		hitlen = read_hits[0] #read_hits[0][-1]

		# divide hit length proportionally among hit taxids
		for i in range(len(all_taxids)):
			this_hitlen = proportions[i] * hitlen
			if args.length_normalization:
				this_hitlen /= tax2info[all_taxids[i]][0]
			if all_taxids[i] in to_add:
				to_add[all_taxids[i]] += this_hitlen
			else:
				to_add[all_taxids[i]] = this_hitlen

	for taxid in to_add:  # add in the multimapped portions
		taxlin = tax2info[taxid][-1].split('|')
		for tax in taxlin:
			if tax != '':
				tax2abs[tax][1] += to_add[taxid]
	return tax2abs


# Renormalize each taxonomic rank so each rank sums to 100% abundance
def rank_renormalize(args, tax2abs, only_strains=False):
	rank_totals = {i:0.0 for i in RANKS}  # current rank abundance sums
	mapped_pct = 100.0
	if not args.no_quantify_unmapped:  # only normalize against the pct of mapped reads
		mapped_pct = 100.0 - (100.0 * tax2abs['Unmapped'][-1])
	for taxid in tax2abs:
		if taxid == 'Unmapped':
			continue
		rank, ab = tax2abs[taxid][1], tax2abs[taxid][-1]
		if only_strains and rank != 'strain':
			continue
		rank_totals[rank] += ab  # add this to the rank sum total

	for taxid in tax2abs:  # here's the normalization
		if taxid == 'Unmapped':
			continue
		rank = tax2abs[taxid][1]
		if only_strains and rank != 'strain':
			continue
		if rank_totals[tax2abs[taxid][1]] != 0.0:
			tax2abs[taxid][-1] /= (rank_totals[tax2abs[taxid][1]] / mapped_pct)
	return tax2abs


# given initital tax2abs, some taxids will be higher than strain level;
#  	we insert "unknown" taxa down to strain level for normalization purposes
def gen_lower_taxa(tax2abs):
	to_add = {}  # lower taxa to add after iteration
	for taxid in tax2abs:
		taxid, rank, taxlin, namelin, ab = tax2abs[taxid]
		if rank == 'strain':  # already lowest level
			continue
		# make a new taxid and name for this strain indicating that it is a
		#  	placeholder for a higher taxon that we cannot further narrow down
		rankpos = RANKS.index(rank)
		lowest_name = namelin.split('|')[rankpos]
		new_name = lowest_name + ' unknown strain'
		new_taxid = taxid+'.0'
		to_add[new_taxid] = [new_taxid, 'strain', taxlin + '|' + new_taxid,
			namelin + '|' + new_name, ab]

	# add in the new taxa
	for taxa in to_add:
		tax2abs[taxa] = to_add[taxa]
	# clean out higher taxa listings -- they will return when we fill the tree
	#tax2abs = {k:v for k,v in tax2abs.items() if v[1] == 'strain'}
	return tax2abs


# Get abundances for all taxids in the tree and put in CAMI format
def format_cami(args, tax2abs):
	# rearrange fields to be in CAMI format
	for taxid in tax2abs:
		old = tax2abs[taxid]
		tax2abs[taxid] = [taxid, old[3], old[5], old[4], old[1]]
	#tax2abs = gen_lower_taxa(tax2abs)  # ensures everything strain level
	# always renormalize strains, to ensure legitimate profile
	#tax2abs = rank_renormalize(args, tax2abs, only_strains=True)
	if not args.no_rank_renormalization:
		tax2abs = rank_renormalize(args, tax2abs)
	return tax2abs


# Performs a Kolmogorov-Smirnov goodness-of-fit test for uniform hit positions
def uniform_startpos_test(acc_hitpos, acclen):
	if len(acc_hitpos) == 0:
		return 'nan', 'nan', 'nan'
	acc_hitpos.sort()
	acc_hitpos = [float(hitpos) / acclen for hitpos, hitlen in acc_hitpos]
	buckets = [0 for i in range(100)]
	for hitpos in acc_hitpos:
		percentile = int(hitpos * 100)
		if percentile >= 100 or percentile <= 0:
			continue
		buckets[percentile] += 1
	test_statistic, pvalue = kstest(acc_hitpos, 'uniform')
	return test_statistic, pvalue, buckets  # acc_hitpos


# Given accession hits, compute the coverage information for each accession
def process_accession_coverages(tax2info, tax2abs, acc2info, acc2hitpos):
	acc2blocks = {}  # collapse hits into blocks of coverage
	acc2basehits = {}  # track total base hits for acc, counting multiple coverage
	acc2unif_test = {}  # track test results for uniformity test (see: uniform_startpos_test)
	for acc in acc2hitpos:
		test_statistic, pvalue, acc_buckets = uniform_startpos_test(acc2hitpos[acc][0], acc2info[acc][0])
		acc2unif_test[acc] = [test_statistic, pvalue, acc_buckets]
		acc2blocks[acc] = [[],[]]
		base_hit_count = 0  # total base hits, including multiple coverage
		for type in range(1):#(2):  # type=0 --> uniq mapped, type=1 --> multi mapped
			hits = acc2hitpos[acc][type]
			if len(hits) == 0:
				continue
			hits.sort(key = lambda x: x[0])  # sort by increasing position
			for hit_start, hit_len in hits:
				hit_start, hit_len = int(hit_start), int(hit_len)
				hit_end = hit_start + hit_len
				base_hit_count += hit_len
				# three cases: create new block, extend it, or don't extend it
				if len(acc2blocks[acc][type]) == 0:
					acc2blocks[acc][type].append([hit_start, hit_end])
				else:
					blocks_end = acc2blocks[acc][type][-1][1]  # end of last block
					if hit_start > blocks_end:
						acc2blocks[acc][type].append([hit_start, hit_end])
					elif hit_end > blocks_end:  # starts in block but extends it
						acc2blocks[acc][type][-1][1] = hit_end
					# last case: hit is subsumed by last block -- do nothing
		acc2basehits[acc] = base_hit_count

	acc2covg = {}
	for acc in acc2blocks:
		acc2covg[acc] = [[],[]]
		for type in range(1):#(2):  # type=0 --> uniq mapped, type=1 --> multi mapped
			if len(acc2blocks[acc][type]) == 0:
				continue
			block_lengths = [block[1] - block[0] for block in acc2blocks[acc][type]]
			covered_bases = sum(block_lengths)
			total_bases = acc2info[acc][0]
			covg_pct = float(covered_bases) / float(total_bases)
			base_hit_count = acc2basehits[acc]
			avg_covg = float(base_hit_count) / float(total_bases)
			tstat, pval, hitpos = acc2unif_test[acc]
			acc2covg[acc][type] = [covg_pct, covered_bases, total_bases,
				base_hit_count, avg_covg, pval, tstat, hitpos]
	return acc2covg


# using the read hit positions in acc2hitpos, compute coverage for accesion for
#  	both uniquely and multi mapped reads. then summarize into coverages for taxids.
def compute_coverages(tax2info, tax2abs, acc2info, acc2hitpos):
	# first use the accession hits to get coverage by accession
	acc2covg = process_accession_coverages(tax2info, tax2abs, acc2info, acc2hitpos)
	# now aggregate these accessions under their TaxID
	tax2covg = {}
	for taxid in tax2info:
		if taxid not in tax2abs:  # a taxid that wasn't mapped to
			continue
		# accs hit, total accs, bases covered, total bases in accs hit & in all accs,
		#   pvalue and test stat. for kstest over whole taxid
		tax2covg[taxid] = [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]]

		taxid_len = tax2info[taxid][0]
		len_processed = 0  # combined lengths of accessions processed so far
		taxid_hitpos = []
		tax2covg[taxid][0][4] = taxid_len
		tax2covg[taxid][1][4] = taxid_len
		for acc in acc2covg:
			if acc2info[acc][1] != taxid:
				continue
			for type in range(1):#(2):  # type=0 --> uniq mapped, type=1 --> multi mapped
				tax2covg[taxid][type][1] += 1
				#tax2covg[taxid][type][4] += acc2info[acc][0]
				if acc in acc2covg and len(acc2covg[acc][type]) != 0:
					tax2covg[taxid][type][0] += 1
					tax2covg[taxid][type][2] += acc2covg[acc][type][1]
					tax2covg[taxid][type][3] += acc2covg[acc][type][2]
					# accumulate accession hit positions for taxid, making sure
					#   that positions are converted to global coordinates
					hitpos = acc2hitpos[acc][type]
					hitpos = [[len_processed + pos[0], pos[1]] for pos in hitpos]
					len_processed += acc2info[acc][0]
					taxid_hitpos.extend(hitpos)
					#tax2covg[taxid][type][5].extend(acc2covg[acc][type][3])
				test_statistic, pvalue, acc_buckets = uniform_startpos_test(taxid_hitpos, taxid_len)
				tax2covg[taxid][type][5:7] = test_statistic, pvalue
	for tax in tax2covg:
		if tax2covg[tax][-1] == 'nan':
			continue
		print(tax, tax2info[tax], '\n', tax2covg[tax][type], '\n')
	return tax2covg, acc2covg


# Mark taxa that are below user-set cutoff thresholds for deletion
def mark_cutoff_taxa(args, tax2info, tax2abs, tax2covg, root = None):
	# test statistics and coverage cutoff values for trimming taxonomic tree up
	#    to a certain level. tstat levels are highest allowed tstat before
	#    trimming that level, wheras coverage levels are lowest allowed.
	tstat_ranks = {'superkingdom': 1.0, 'phylum': 0.8, 'class': 0.6,
		'order': 0.4, 'family': 0.3, 'genus': 0.2, 'species': 0.1, 'strain': 0.01}
	covg_ranks = {'superkingdom': 0.001, 'phylum': 0.001, 'class': 0.001,
		'order': 0.005, 'family': 0.01, 'genus': 0.05, 'species': 0.1, 'strain': 0.5}

	del_list, keep_list = {}, {}  # taxa to delete or keep from the tree
	for taxid in tax2abs:
		if taxid == 'Unmapped' or 'Unmapped' in taxid:
			continue
		taxinfo = tax2abs[taxid]
		num_reads, taxlin = taxinfo[0], taxinfo[-1]
		taxlin_splits = taxlin.split('|')
		if num_reads < args.read_cutoff:
			del_list[taxid] = True
			continue
		if taxid in tax2covg:
			uniq_stats = tax2covg[taxid][0]
			t_stat, pval = float(uniq_stats[5]), float(uniq_stats[6])
			bases_covd, total_bases = uniq_stats[2], uniq_stats[3]
			if total_bases == 0:
				uniq_covg_pct = 0.0
			else:
				uniq_covg_pct = float(bases_covd) / float(total_bases)

			# here we check how far up the branch to trim -- see description above
			t_stat_trimlevel, covg_trimlevel = 1, 1
			for i in range(1, len(RANKS) + 1):
				rank_to_check = RANKS[-i]  # rank to check whether to trim
				if t_stat > tstat_ranks[rank_to_check]:
					t_stat_trimlevel = i
				if uniq_covg_pct < covg_ranks[rank_to_check]:
					covg_trimlevel = i
			# note: we take the stricter of the two trims
			overall_trimlevel = max(t_stat_trimlevel, covg_trimlevel)
			print(taxid, overall_trimlevel, t_stat, uniq_covg_pct)

			for i in range(1, len(RANKS) + 1):
				cur_tax = taxlin_splits[-i]
				if i <= overall_trimlevel and cur_tax not in keep_list:
					del_list[cur_tax] = True
				elif i > overall_trimlevel:
					keep_list[cur_tax] = True
					if cur_tax in del_list:
						del del_list[cur_tax]
			#if uniq_covg_pct < args.uniq_covg_cutoff:
			#	del_list[taxid] = True
	return del_list


# ensure the children of pruned parents are also pruned, as well as childless parents
def ensure_consistency(tax2info, tax2abs, del_list):
	# ensure the children of pruned parents/ancestors are also pruned
	for taxid in tax2abs:
		taxlin = tax2abs[taxid][-1]
		taxlin_splits = taxlin.split('|')
		if len(taxlin_splits) == 1:  # superkingdom level should have no parent
			continue
		for i in range(2, len(taxlin_splits) + 1):
			ancestor = taxlin.split('|')[-i]
			if ancestor == '':
				continue
			if ancestor in del_list:
				del_list[taxid] = True
				break
	# prune nodes with no descendents
	'''for i in range(len(RANKS)):  # do for each taxonomic rank
		for taxid in tax2abs:
			if tax2abs[taxid][3] != RANKS[len(RANKS) - i -1]:
				continue
			if taxid not in tax2info:  # not leaf node
				taxlin = tax2abs[taxid][-1]
				descendents = [node for node in tax2abs if taxlin + '|' in tax2abs[node][-1]]
				descendents = [desc for desc in descendents if desc not in del_list]
				if len(descendents) == 0:
					del_list[taxid] = True'''
	return del_list


# After marking all deletable taxa, see if parents of pruned children should survive
#  when considering reads multimapped for children but unique mapped for parents
def rescue_higher_taxa(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos, del_list):
	for taxid in tax2abs:
		if taxid in tax2info or taxid not in del_list:
			continue  # skip leaf nodes or those not facing deletion
		# hits that are unique to this taxa but multimapped in lower taxa
		multi_hits_in_lower = tax2abs[taxid][2]
		new_acc2hitpos = {}  # count those reads as unique in leaf nodes
		for acc in multi_hits_in_lower:
			new_acc2hitpos[acc] = [[],[]]
			new_acc2hitpos[acc][0] = acc2hitpos[acc][0]
			# indices of reads multimapped at lowest level that will count as unique now
			indices_to_add = [acc2hitpos[acc][1][i] for i in multi_hits_in_lower[acc]
				if len(acc2hitpos[acc][1]) != 0]
			new_acc2hitpos[acc][0].extend(indices_to_add)

		# generate this node's subtree, recompute coverages, and see if it would be pruned
		this_tax2abs = {}
		taxlin = tax2abs[taxid][-1]
		this_tax2abs = {node: tax2abs[node] for node in tax2abs
			if taxlin + '|' in tax2abs[node][-1]}
		this_tax2abs[taxid] = tax2abs[taxid]  # put current node in there
		new_tax2covg = compute_coverages(tax2info, this_tax2abs, acc2info, new_acc2hitpos)
		new_del_list = mark_cutoff_taxa(args, tax2info, this_tax2abs, new_tax2covg, root = taxid)
		new_del_list = ensure_consistency(tax2info, this_tax2abs, new_del_list)
		if taxid not in new_del_list:
			del del_list[taxid]
	return del_list


# Apply read cutoff and min_abundance thresholds, but latter relative to parent taxa
def prune_tree(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos):
	del_list = mark_cutoff_taxa(args, tax2info, tax2abs, tax2covg)
	del_list = ensure_consistency(tax2info, tax2abs, del_list)
	#del_list = rescue_higher_taxa(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos, del_list)

	# delete taxids with insufficient evidence, keeping their reads in an
	#    "Unmapped" entry for abundance estimation purposes
	for taxid in del_list:
		#rank = tax2abs[taxid][-3]
		#unm_taxid = rank + '_Unmapped'
		#if unm_taxid not in tax2abs:
		#	tax2abs[unm_taxid] = tax2abs[taxid]
		#else:
		#	tax2abs[unm_taxid][0] += tax2abs[taxid][0]
		#	tax2abs[unm_taxid][1] += tax2abs[taxid][1]
		del tax2abs[taxid]
	return tax2abs


# Processes uniquely-mapped reads, then estimates abundances using
#  	uniquely-mapped abundances and multimapped read information
def compute_abundances(args, infile, acc2info, tax2info, acc2uniqblocks):
	echo('Reading input file ' + infile, args.verbose)
	# run mapping and process to get uniq map abundances & multimapped reads
	tax2abs, multimapped, acc2hitpos = map_and_process(args, infile,
		acc2info, tax2info)

	echo('Breaking accessions into unique blocks...', args.verbose)
	acc2hitpos, acc2info, tax2info, tax2abs = break_accs_into_uniq_blocks(
		acc2hitpos, acc2info, tax2info, tax2abs, acc2uniqblocks)

	echo('Computing coverage information...', args.verbose)
	tax2covg, acc2covg = compute_coverages(tax2info, tax2abs, acc2info, acc2hitpos)
	# prune tree based on user cutoff settings
	tax2abs = prune_tree(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos)

	echo('Assigning multimapped reads...', args.verbose)
	tax2abs = resolve_multi_prop(args, tax2abs, multimapped, tax2info)
	results = format_cami(args, tax2abs)
	echo('Done computing abundances for input file ' + infile, args.verbose)
	return results


# Combines and averages results across all input files,
#  	and packs information into easy-to-write form
def gather_results(args, acc2info, acc2uniqblocks, tax2info):
	results = {}
	for infile in args.infiles:
		echo('Computing abundances for input file ' + infile + '...', args.verbose)
		file_res = compute_abundances(args, infile, acc2info, tax2info, acc2uniqblocks)
		for taxid in file_res:
			if taxid not in results:
				results[taxid] = file_res[taxid]
			else:
				results[taxid][-1] += file_res[taxid][-1]

	echo('Compiling and writing results...', args.verbose)
	rank_results = {i:[] for i in range(len(RANKS))}
	for taxid in results:
		if 'Unmapped' in taxid:
			continue
		ab = results[taxid][-1]  # avg over input files, truncate to 4 digits
		if results[taxid][-1] < 0.00001:
			results[taxid][-1] = '0.00001'
		else:
			results[taxid][-1] = float('%.5f' % results[taxid][-1])
		#results[taxid][-1] = float(math.trunc(
		#	(ab / len(args.infiles)) * (10**4))) / (10**4)
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
			outfile.write('@SampleID:' + ','.join(args.infiles) + '\n')
		else:
			outfile.write('@SampleID:' + args.sampleID + '\n')
		outfile.write('@Version:MiCoP2-v0.1\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t' +
			'PERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

		for i in range(len(RANKS)):
			if args.no_strain_level and i == len(RANKS)-1:
				continue  # skip the strain level if user wants to
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


def map_main(args = None):
	if args == None:
		args = profile_parseargs()
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		print('Error: --pct_id must be between 0.0 and 1.0, inclusive.')
		sys.exit()
	if args.db == 'NONE' and not args.infiles[0].endswith('sam'):
		print('Error: --db must be specified unless .sam files are provided.')
		sys.exit()
	if args.dbinfo == 'AUTO':
		args.dbinfo = __location__ + 'data/subset_db_info.txt'
	open(args.output, 'w').close()  # test to see if writeable

	# maps NCBI accession to length, taxid, name lineage, taxid lineage
	echo('Reading dbinfo and uniq_blocks files...', args.verbose)
	acc2info, tax2info = get_acc2info(args)
	acc2uniqblocks = read_unique_blocks(args.uniq_blocks_file, acc2info)
	# gathers results for all infiles, combines, and organizes into tax levels
	rank_results = gather_results(args, acc2info, acc2uniqblocks, tax2info)
	write_results(args, rank_results)


if __name__ == '__main__':
	args = profile_parseargs()
	map_main(args)
#
