import argparse, math, os, subprocess, sys, time


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


def profile_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Compute abundance estimations for species in a sample.')
	parser.add_argument('infiles', nargs='+',
		help='sam or reads file(s) (space-delimited if multiple). Required.')
	parser.add_argument('--db', default='NONE',
		help='Path to database from select_db.py. Required if read files given.')
	parser.add_argument('--dbinfo', default='AUTO',
		help = 'Location of db_info file. Default: data/subset_db_info.txt')
	parser.add_argument('--min_abundance', type=float, default=10**-4,
		help='Minimum abundance for a taxa to be included in the results.')
	parser.add_argument('--length_normalization', action='store_true',
		help='Normalize abundances by genome length.')
	parser.add_argument('--no_rank_renormalization', action='store_true',
		help='Do not renormalize abundances to 100 percent at each rank,\
				for instance if an organism has a species but not genus label.')
	parser.add_argument('--output', default='abundances.tsv',
		help='Output abundances file. Default: abundances.txt')
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
	parser.add_argument('--uniq_covg_cutoff', default = 0.2, type = float,
		help='Coverage pct. by uniquely-mapped reads needed to count as present.')
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
def get_lowest_common_ancestor(read_hits, taxid2info):
	lca = ''
	all_taxids = [hit[7] for hit in read_hits]
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


# Given hit information and LCA, update abundance info for relevant taxa
def process_hit(hit, taxid, lca_index, avg_hitlen, acc2hitpos, tax2abs, taxid2info):
	taxids_processed = {}
	# record position and hitlen for lowest level taxid hit
	acc, pos, pct_id, hitlen = hit[2], hit[3], hit[5][0], hit[5][1]
	if acc not in acc2hitpos:
		acc2hitpos[acc] = [[],[]]
	if taxid != 'multi':
		acc2hitpos[acc][0].append([pos, hitlen])
	else:
		acc2hitpos[acc][1].append([pos, hitlen])

	taxlin = taxid2info[hit[7]][-1].split('|')
	namelin = taxid2info[hit[7]][-2].split('|')
	for i in range(len(taxlin)):
		this_taxid = taxlin[i]
		if this_taxid == '' or this_taxid in taxids_processed:
			continue
		taxids_processed[this_taxid] = ''  # avoid re-processing this taxid
		is_unique = (i <= lca_index)  # unique hit if LCA or an ancestor
		if this_taxid in tax2abs:
			if this_taxid in taxid2info:  # lowest level taxid
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
			if this_taxid in taxid2info:  # lowest level taxid
				# 3rd field is accession length -- for covg purposes
				if is_unique:
					tax2abs[this_taxid] = [1, avg_hitlen] + taxid2info[this_taxid]
				#else:
				#	tax2abs[this_taxid] = [0, 0] + taxid2info[this_taxid]
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
def update_global_hits(args, hits, taxid, global_results, taxid2info, acc2info):
	tax2abs, multimapped, acc2hitpos = global_results
	if taxid == 'Ambiguous':  # ambiguous mapping (see intersect_read_hits)
		if not args.no_quantify_unmapped:
			if 'Unmapped' in tax2abs:
				tax2abs['Unmapped'][0] += 1.0  # reads hit
			else:  # also store lineage for taxid
				tax2abs['Unmapped'] = ([1.0, 0.0] +
					taxid2info['Unmapped'])
		return tax2abs, multimapped, acc2hitpos
	if taxid != 'multi' and args.length_normalization:
		for i in range(len(hits)):
			hits[i][5][1] /= taxid2info[hits[i][7]][0] # normalize by genome length
	lca = get_lowest_common_ancestor(hits, taxid2info)
	if lca == '':
		lca_index = -1
	else:
		lca_index = taxid2info[hits[0][7]][-1].split('|').index(lca)
	avg_hitlen = float(sum([hit[5][1] for hit in hits])) / len(hits)
	for hit in hits:
		tax2abs, taxids_processed, acc2hitpos = process_hit(hit, taxid, lca_index,
			avg_hitlen, acc2hitpos, tax2abs, taxid2info)

	if taxid == 'multi':  # store taxids hit and length (in bases) of hits
		#hits = [[hit[7], len(hit[9])] for hit in hits]
		hits = [avg_hitlen] + [hit[7] for hit in hits]
		multimapped.append(hits)
	return tax2abs, multimapped, acc2hitpos


# Runs minimap2, processes output, and returns abundance information for taxids,
#  and multimapping and coverage information for post-processing.
# This is the general iterator function; it calls subroutines to process
#  individual lines of input
def map_and_process(args, infile, acc2info, taxid2info):
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
					tax2abs['Unmapped'] = [1.0, 0.0] + taxid2info['Unmapped']
			continue
		splits = line.strip().split()
		pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
		if is_bad:  # unmapped or cigar string unavailable
			continue
		splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
		if '|' in line:
			acc = splits[2].split('|')[-2]
			splits[2] = acc
		else:
			acc = splits[2]
		splits[7] = acc2info[acc][1]  # replace PNEXT field (unused) w/ taxid

		read = splits[0]
		if read != prev_read:
			total_reads += 1
			if args.verbose and total_reads % 100000 == 0:
				echo('Done processing ' + str(total_reads) + ' reads.')
			# get uniq hit taxid or multimapped hits intersect
			read_hits, pair1maps, pair2maps = filter_read_hits(
				args, read_hits, pair1maps, pair2maps)
			intersect_hits, taxid = intersect_read_hits(
				args, read_hits, pair1, pair2, pair1maps, pair2maps)
			# update accumulated global results using the read hits information
			global_results = [tax2abs, multimapped, acc2hitpos]
			tax2abs, multimapped, acc2hitpos = update_global_hits(
				args, intersect_hits, taxid, global_results, taxid2info, acc2info)
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


# Divide abundances of multimapped reads according to uniquely mapped reads
#  	portion for each of the hit organisms
def resolve_multi_prop(args, tax2abs, multimapped, taxid2info):
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
				this_hitlen /= taxid2info[all_taxids[i]][0]
			if all_taxids[i] in to_add:
				to_add[all_taxids[i]] += this_hitlen
			else:
				to_add[all_taxids[i]] = this_hitlen

	for taxid in to_add:  # add in the multimapped portions
		taxlin = taxid2info[taxid][-1].split('|')
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


# Given accession hits, compute the coverage information for each accession
def process_accession_coverages(taxid2info, tax2abs, acc2info, acc2hitpos):
	acc2blocks = {}  # collapse hits into blocks of coverage
	for acc in acc2hitpos:
		acc2blocks[acc] = [[],[]]
		for type in range(1):#(2):  # type=0 --> uniq mapped, type=1 --> multi mapped
			hits = acc2hitpos[acc][type]
			if len(hits) == 0:
				continue
			hits.sort(key = lambda x: x[0])  # sort by increasing position
			for hit_start, hit_len in hits:
				hit_start, hit_len = int(hit_start), int(hit_len)
				hit_end = hit_start + hit_len
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

	acc2covg = {}
	for acc in acc2blocks:
		acc2covg[acc] = [[],[]]
		for type in range(1):#(2):  # type=0 --> uniq mapped, type=1 --> multi mapped
			if len(acc2blocks[acc][type]) == 0:
				continue
			block_lengths = [block[1] - block[0] for block in acc2blocks[acc][type]]
			covered_bases = sum(block_lengths)
			total_bases = acc2info[acc][0]
			coverage = float(covered_bases) / float(total_bases)
			acc2covg[acc][type] = [coverage, covered_bases, total_bases]#, block_lengths]
	return acc2covg


# using the read hit positions in acc2hitpos, compute coverage for accesion for
#  	both uniquely and multi mapped reads. then summarize into coverages for taxids.
def compute_coverages(taxid2info, tax2abs, acc2info, acc2hitpos):
	# first use the accession hits to get coverage by accession
	acc2covg = process_accession_coverages(taxid2info, tax2abs, acc2info, acc2hitpos)
	# now aggregate these accessions under their TaxID
	tax2covg = {}
	for taxid in taxid2info:
		if taxid not in tax2abs:  # a taxid that wasn't mapped to
			continue
		# accs hit, total accs, bases covered, total bases in accs hit & in all accs, block lengths
		tax2covg[taxid] = [[0, 0, 0, 0, 0, []], [0, 0, 0, 0, 0, []]]

		taxid_len = taxid2info[taxid][0]
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
					#tax2covg[taxid][type][5].extend(acc2covg[acc][type][3])
	return tax2covg


# Mark taxa that are below user-set cutoff thresholds for deletion
def mark_cutoff_taxa(args, tax2info, tax2abs, tax2covg, root = None):
	del_list = {}
	for taxid in tax2abs:
		if taxid == 'Unmapped':
			continue
		taxinfo = tax2abs[taxid]
		num_reads, taxlin = taxinfo[0], taxinfo[-1]
		if num_reads < args.read_cutoff:
			del_list[taxid] = True
			continue
		if taxid in tax2covg:
			uniq_stats = tax2covg[taxid][0]
			bases_covd, total_bases = uniq_stats[2], uniq_stats[3]
			if total_bases == 0:
				uniq_covg_pct = 0.0
			else:
				uniq_covg_pct = float(bases_covd) / float(total_bases)
			if uniq_covg_pct < args.uniq_covg_cutoff:
				del_list[taxid] = True

		taxlin_splits = taxlin.split('|')
		if len(taxlin_splits) == 1 or taxid == root:  # no parent
			continue
		for i in range(2, len(taxlin_splits) + 1):
			parent = taxlin.split('|')[-i]
			if parent != '':
				break
		parent_reads = tax2abs[parent][0]
		if float(num_reads) / float(parent_reads) < args.min_abundance:
			del_list[taxid] = True
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
	for i in range(len(RANKS)):  # do for each taxonomic rank
		for taxid in tax2abs:
			if tax2abs[taxid][3] != RANKS[len(RANKS) - i -1]:
				continue
			if taxid not in tax2info:  # not leaf node
				taxlin = tax2abs[taxid][-1]
				descendents = [node for node in tax2abs if taxlin + '|' in tax2abs[node][-1]]
				descendents = [desc for desc in descendents if desc not in del_list]
				if len(descendents) == 0:
					del_list[taxid] = True
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
	del_list = rescue_higher_taxa(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos, del_list)
	for taxid in del_list:
		del tax2abs[taxid]
	return tax2abs


# Processes uniquely-mapped reads, then estimates abundances using
#  	uniquely-mapped abundances and multimapped read information
def compute_abundances(args, infile, acc2info, tax2info):
	if args.verbose:
		echo('Reading input file ' + infile)
	# run mapping and process to get uniq map abundances & multimapped reads
	tax2abs, multimapped, acc2hitpos = map_and_process(args, infile,
		acc2info, tax2info)
	if args.verbose:
		echo('Computing coverage information...')
	tax2covg = compute_coverages(tax2info, tax2abs, acc2info, acc2hitpos)
	# prune tree based on user cutoff settings
	tax2abs = prune_tree(args, tax2info, tax2abs, tax2covg, acc2info, acc2hitpos)

	if args.verbose:
		echo('Assigning multimapped reads...')
	tax2abs = resolve_multi_prop(args, tax2abs, multimapped, tax2info)
	results = format_cami(args, tax2abs)
	if args.verbose:
		echo('Multimapped reads assigned.')
		echo('Done computing abundances for input file ' + infile)
	return results


# Combines and averages results across all input files,
#  	and packs information into easy-to-write form
def gather_results(args, acc2info, taxid2info):
	results = {}
	for infile in args.infiles:
		if args.verbose:
			echo('Computing abundances for input file ' + infile + '...')
		file_res = compute_abundances(args, infile, acc2info, taxid2info)
		for taxid in file_res:
			if taxid not in results:
				results[taxid] = file_res[taxid]
			else:
				results[taxid][-1] += file_res[taxid][-1]

	if args.verbose:
		echo('Compiling and writing results...')
	rank_results = {i:[] for i in range(len(RANKS))}
	for taxid in results:
		ab = results[taxid][-1]  # avg over input files, truncate to 4 digits
		if results[taxid][-1] < 0.00001:
			results[taxid][-1] = '0.00001'
		else:
			results[taxid][-1] = str(float('%.5f' % results[taxid][-1]))
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
				if line[4] < 0.0001:  # args.min_abundance:
					continue
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
	acc2info, taxid2info = get_acc2info(args)
	# gathers results for all infiles, combines, and organizes into tax levels
	rank_results = gather_results(args, acc2info, taxid2info)
	write_results(args, rank_results)


if __name__ == '__main__':
	args = profile_parseargs()
	map_main(args)
#
