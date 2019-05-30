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


def profile_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Compute abundance estimations for species in a sample.')
	parser.add_argument('infiles', nargs='+',
		help='sam or reads file(s) (space-delimited if multiple). Required.')
	parser.add_argument('--assignment', choices=['em', 'proportional', 'none'],
		default='proportional', help='Method for assignming multimapped reads.')
	parser.add_argument('--db', default='NONE',
		help='Path to database from select_db.py. Required if read files given.')
	parser.add_argument('--dbinfo', default='AUTO',
		help = 'Location of db_info file. Default: data/subset_db_info.txt')
	parser.add_argument('--min_abundance', type=float, default=10**-4,
		help='Minimum abundance for a taxa to be included in the results.')
	parser.add_argument('--min_map', type=int, default=-1,
		help='Minimum bases mapped to count a hit.')
	parser.add_argument('--max_ed', type=int, default=999999999,
		help='Maximum edit distance from a reference to count a hit.')
	parser.add_argument('--no_len_normalization', action='store_true',
		help='Do not normalize abundances by genome length.')
	parser.add_argument('--no_rank_renormalization', action='store_true',
		help='Do not renormalize abundances to 100 percent at each rank,\
				for instance if an organism has a species but not genus label.')
	parser.add_argument('--output', default='abundances.tsv',
		help='Output abundances file. Default: abundances.txt')
	parser.add_argument('--pct_id', type=float, default=-1,
		help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--quantify_unmapped', action='store_true',
		help='Factor in unmapped reads in abundance estimation.')
	parser.add_argument('--read_cutoff', type=int, default=-1,
		help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--strain_level', action='store_true',
		help='Write output at the strain level as well.')
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


# Determine whether to count a mapping hit based on user specifications
# Return True if hit should be filtered out, False otherwise
def filter_line(args, splits):
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
	if matched_len < args.min_map:
		return True
	edit_distance = int(splits[11][5:])
	if edit_distance > args.max_ed:
		return True
	#if float(mapped) / float(mapped + edit_distance) < args.pct_id:
	if float(matched_len) / float(total_len) < args.pct_id:
		return True
	return False  # if read passes quality checks, don't filter


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
		# remove user-filtered & chimeric reads, update # of pair end map counts
		if filter_line(args, read_hits[hit]) or read_hits[hit][1][2]:
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
			return [], read_hits[0][2], readquals, hitlen

		intersect, intersect_hits = intersect_read_hits(
			read_hits, pair1maps, pair2maps)
		if len(intersect) == 0:  # one end unmapped, other multimapped
			return [], 'Ambiguous', '', -1  #we consider this case too ambiguous
		elif len(intersect) == 1:  # read pairs only agree on one taxid
			return [], read_hits[0][2], readquals, hitlen  # considered uniq map
		else:  # both ends multimapped, handle later
			return intersect_hits, '', readquals, hitlen

	else:  # single end
		if pair1maps > 1:  # multimapped
			return read_hits, '', readquals, hitlen  # multimapped
		else:
			return False, read_hits[0][2], readquals, hitlen


# Computes the likelihood estimate for a read mapping using base quality scores
#  	and CIGAR string
def single_read_likelihood(probs, cigar):
	# curval = current cigar letter amount, cur_ind = base probs. index
	base_likelihoods, curval, cur_ind = [], 0, 0
	for ch in cigar:
		if ch.isdigit():
			curval = (curval * 10) + int(ch)
		else:
			if ch != 'D':  # deletions do not consume base quality scores
				if ch == 'M' or ch == '=':  # base likelihoods for matched bases
					base_likelihoods.extend(
						[1-probs[i] for i in range(cur_ind, cur_ind + curval)])
				else:  # base likelihoods for mismatched base
					base_likelihoods.extend(
						[probs[i]/3.0 for i in range(cur_ind, cur_ind+curval)])
				cur_ind = cur_ind + curval
			curval = 0
	return reduce(mul, base_likelihoods, 1)  # product of base likelihoods


# computes likelihood of read being assigned to each hit taxid based on CIGAR
#  	string and base quality scores
def compute_read_likelihoods(read_hits):
	base_quals, hitlen = read_hits[0][10], len(read_hits[0][10])
	hit_dict = {}  # used to group reads by reference hit, needed for paired end
	for hit in read_hits:
		if hit[2] in hit_dict:
			hit_dict[hit[2]].append(hit)
		else:
			hit_dict[hit[2]] = [hit]

	likelihoods = {}
	for taxid in hit_dict:
		cigars = hit_dict[taxid][0][5]
		if len(hit_dict[taxid]) == 2:  # paired end
			cigars += hit_dict[taxid][1][5]
		else:
			# eliminates unfair advantages for one end mapped reads caused by
			#  	having less probabilities to multiply over than both ends mapped
			cigars += hit_dict[taxid][0][5]
		# compute prob. of each base being wrong, then this read's likelihood
		prob_bases_wrong = [10 ** (-(ord(ch) - 33) / 10) for ch in base_quals]
		read_likelihood = single_read_likelihood(prob_bases_wrong, cigars)
		if read_likelihood  < 10 ** -300:  # minimum likelihood filter
			likelihoods[taxid] = read_likelihood
	return [likelihoods, hitlen]


# Remove hits to taxids not in taxids2abs (no unique mappings to that taxid)
def preprocess_multimapped(args, total_bases, mmap_bases, multimapped, taxids2abs):
	for i in range(len(multimapped)):
		if args.assignment == 'proportional':
			hitlen = multimapped[i][0][-1]
			multimapped[i] = [hit for hit in multimapped[i]
				if hit[0] in taxids2abs]
			if len(multimapped[i]) > 0 and len(multimapped[i][0]) == 1:
				multimapped[i][0].append(hitlen)  # ensure hitlen is kept
		else:
			multimapped[i][0] = {hit: val
				for hit,val in multimapped[i][0].items() if hit in taxids2abs}
		if len(multimapped[i]) > 0:
			total_bases += mmap_bases[i]
	multimapped = [read for read in multimapped if len(read) > 0]
	return total_bases, multimapped


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


# Runs minimap2, processes output, and returns a dict of taxids mapped to number
#  	of uniquely mapped reads and bases for that taxid,
#  	and a list of multimapped reads.
def map_and_process(args, infile, acc2info, taxid2info):
	taxids2abs, multimapped = {}, []  # taxids to abundances, multimapped reads
	mmap_bases = []  # tracks number of bases for each multimapped read
	prev_read, read_hits = '', [] # read tracker, all hits for read (full lines)
	pair1maps, pair2maps = 0, 0  # reads mapped to each pair (single = pair1)
	total_bases, total_reads = 0, 0

	samfile = False  # whether reading from existing sam file
	if infile.endswith('sam'):  # input stream from sam file
		samfile = True
		instream = open(infile, 'r')
	else:  # run minimap2 and stream its output as input
		mapper = subprocess.Popen(['./minimap2/minimap2', '-ax', 'sr',
			'-t', '4', '-2', '-n' '1', '--secondary=yes',
			args.db, infile], stdout=subprocess.PIPE, bufsize=1)
		instream = iter(mapper.stdout.readline, "")

	for line in instream:
		if not samfile:
			line = line.decode('utf-8')
			if not line:  # process finished
				break
		if line.startswith('@'):
			if args.quantify_unmapped:
				if 'Unmapped' in taxids2abs:
					taxids2abs['Unmapped'][0] += 1.0  # reads hit
				else:  # also store lineage for taxid
					taxids2abs['Unmapped'] = [1.0, 0.0] + taxid2info['Unmapped']
			continue
		splits = line.strip().split()
		pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
		if is_bad:  # unmapped or cigar string unavailable
			continue
		splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
		# here we change accession to taxid since we want to
		#  	profile by organisms, not accessions
		#splits[2] = acc2info[splits[2]][1]
		acc = splits[2].split('|')[-2]
		splits[2] = acc2info[acc][1]

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
				if args.quantify_unmapped:
					if 'Unmapped' in taxids2abs:
						taxids2abs['Unmapped'][0] += 1.0  # reads hit
					else:  # also store lineage for taxid
						taxids2abs['Unmapped'] = ([1.0, 0.0] +
							taxid2info['Unmapped'])
				continue
			if taxid != '' and not (args.no_len_normalization
					or args.assignment == 'em'):
				hitlen /= taxid2info[taxid][0]  # normalize by genome length
			if intersect_hits == []:  # unique hit
				total_bases += hitlen
				taxlin = taxid2info[taxid][-1].split('|')
				namelin = taxid2info[taxid][-2].split('|')
				#taxid2info[taxid] = [acclen, rank, namelin, taxlin]
				for i in range(len(taxlin)):
					this_taxid = taxlin[i]
					if this_taxid == '':
						continue
					if this_taxid in taxids2abs:
						taxids2abs[this_taxid][0] += 1
						taxids2abs[this_taxid][1] += hitlen
					else:
						this_rank = RANKS[i]
						this_taxlin = '|'.join(taxlin[:i+1])
						this_namelin = '|'.join(namelin[:i+1])
						taxids2abs[this_taxid] = [1, hitlen] + [0,
							this_rank, this_namelin, this_taxlin]
				'''if taxid in taxids2abs:
					taxids2abs[taxid][0] += 1  # reads hit
					taxids2abs[taxid][1] += hitlen  # bases hit
				else:  # also store lineage for taxid
					taxids2abs[taxid] = [1, hitlen] + taxid2info[taxid]'''
			else:  # multimapped hit
				lca = get_lowest_common_ancestor(intersect_hits, taxid2info)
				if lca != '':
					# counts as unique hit for the LCA and higher ranks
					taxinfo = taxid2info[intersect_hits[0][2]]
					full_taxlin = taxinfo[-1].split('|')
					full_namelin = taxinfo[-2].split('|')
					lca_index = full_taxlin.index(lca)
					taxlin = full_taxlin[:lca_index + 1]
					namelin = full_taxlin[:lca_index + 1]
					#taxlin = taxid2info[lca][-1].split('|')
					#namelin = taxid2info[lca][-2].split('|')
					#taxid2info[taxid] = [acclen, rank, namelin, taxlin]
					for i in range(len(taxlin)):
						this_taxid = taxlin[i]
						if this_taxid == '':
							continue
						if this_taxid in taxids2abs:
							taxids2abs[this_taxid][0] += 1
							taxids2abs[this_taxid][1] += hitlen
						else:
							this_rank = RANKS[i]
							this_taxlin = '|'.join(taxlin[:i+1])
							this_namelin = '|'.join(namelin[:i+1])
							taxids2abs[this_taxid] = [1, hitlen] + [0,
								this_rank, this_namelin, this_taxlin]

				if args.assignment == 'proportional':
					# store taxids hit and length (in bases) of hits
					#intersect_hits = [[hit[2], len(hit[9])]
					#	for hit in intersect_hits]
					intersect_hits = [[hit[2]] for hit in intersect_hits]
					intersect_hits[0].append(len(readquals))  # total hit length
				elif args.assignment == 'em':  # compute assingment likelihoods
					intersect_hits[0][10] = readquals  # ensure qual scores
					intersect_hits[0].append(pair1 or pair2)  # paired or not
					intersect_hits = compute_read_likelihoods(intersect_hits)
				# for em, ensure min. likelihood filter didn't filter all hits
				if intersect_hits[0] != {} and args.assignment != 'none':
					multimapped.append(intersect_hits)
					mmap_bases.append(len(readquals))
		#else:
		pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
		pair2maps += pair2  # unchanged if pair2 false
		read_hits.append(splits)

	if samfile:
		instream.close()
	else:
		mapper.stdout.close()
		mapper.wait()
	if len(multimapped) > 0:
		unique_bases = total_bases
		total_bases, multimapped = preprocess_multimapped(args, total_bases,
			mmap_bases, multimapped, taxids2abs)
		mmap_bases = total_bases - unique_bases  # number of multimapped bases
	else:
		mmap_bases = 0
	# store percentage of unmapped reads
	if args.quantify_unmapped:
		taxids2abs['Unmapped'][1] = taxids2abs['Unmapped'][0] / float(total_reads)
	return taxids2abs, multimapped, float(total_bases), float(mmap_bases)


# Initial EM abundances estimate: similar to Karp paper, except proportion of
#  	uniquely-mapped reads is changed to (possibly length-normalized) proportion
#  	of uniquely-mapped bases, allowing variable read length to factor in
def initital_estimate(taxids2abs, multimapped, total_bases, mmap_bases):
	# calculate share of total bases uniquely mapped to each taxid
	taxids2abs = {k: v/total_bases for k,v in taxids2abs.items()}
	uniq_proportions = {k:v for k,v in taxids2abs.items()}

	# now add in initial estimate for multimapped reads: evenly-distributed
	all_mmap_taxids = []  # set of multimapped taxids with >= 1 uniq mapped read
	for read in multimapped:
		all_mmap_taxids.extend([hit for hit in read[0] if hit in taxids2abs])
	all_mmap_taxids = set(all_mmap_taxids)
	even_dist = mmap_bases / total_bases / float(len(all_mmap_taxids))
	for taxid in taxids2abs:
		taxids2abs[taxid] += even_dist
	return uniq_proportions, taxids2abs


# Based on magnitude of changes since last EM estimate, decide whether to stop
def end_condition(changes, prev_magnitude):
	ab_magnitude = sum([abs(v) for k,v in changes.items()])
	if ab_magnitude < (10 ** -3):
		return True, ab_magnitude
	diff_pct = abs(prev_magnitude - ab_magnitude) / prev_magnitude
	if diff_pct < 0.1:
		return True, ab_magnitude
	return False, ab_magnitude


# Performs an updated abundance estimation for an interation of the EM algorithm
def em_update(args, total_bases, uniq_prop, taxids2abs,taxid2info, multimapped):
	init_abs = {k:v for k,v in taxids2abs.items()}
	# get read assignment likelihoods for each read
	likelihood_sums = {}
	for read in multimapped:
		read_likelihoods, hitlen = read  # unpack read info
		# base_share is the share of total bases reflected in this read
		base_share = float(hitlen) / float(total_bases)
		# calculate read likelihoods times abundances ("ab_likelihoods") and
		#  	the sum of the ab_likelihoods, and the ratio of each to the total
		ab_likelihoods = {k: v * taxids2abs[k] for k,v in read_likelihoods.items()
			if k in taxids2abs}
		total_likelihood = sum([v for k,v in ab_likelihoods.items()])
		if total_likelihood == 0.0:
			continue
		read_ratios = {k: v / total_likelihood for k,v in ab_likelihoods.items()}
		# sum all likelihood ratios for a taxid (combines multiple hits)
		for taxid in read_ratios:
			if taxid not in likelihood_sums:
				likelihood_sums[taxid] = (read_ratios[taxid] * base_share)
			else:  # add in additional hits for the same taxid
				likelihood_sums[taxid] += (read_ratios[taxid] * base_share)

	# update taxids2abs
	for taxid in likelihood_sums:#taxids2abs:
		taxids2abs[taxid] = uniq_prop[taxid] + likelihood_sums[taxid]
	if not args.no_len_normalization:  # normalize by genome len. if user wants
		taxids2abs = {k: v / float(taxid2info[k][0])
			for k,v in taxids2abs.items()}
		all_abs = sum([v for k,v in taxids2abs.items()])
		taxids2abs = {k: v/all_abs for k,v in taxids2abs.items()}  # renorm
	changes = {k:(taxids2abs[k]-v) for k,v in init_abs.items()}
	return changes, taxids2abs


# Use EM algorithm to revise taxid abundances using multimapped information,
#  	without explicitly assigning reads. Recommended method in MiCoP2.
def resolve_multi_em(args, total_bases, mmap_bases, taxids2abs, multimapped, taxid2info):
	# we temporarily just keep abundance to simplify EM operations
	only_abs = {k: v[1] for k,v in taxids2abs.items() if k != 'Unmapped'}
	if len(multimapped) == 0:
		return taxids2abs
	uniq_prop, only_abs = initital_estimate(
		only_abs, multimapped, total_bases, mmap_bases)

	# we have a dict showing changes in ab. estimates for each step;
	#  	 if this changes dict has small enough magnitude, end EM updates
	changes = {k:(only_abs[k]-v) for k,v in uniq_prop.items()}
	em_iter, end_iters = 0, 0
	end, prev_magnitude = end_condition(changes, 1000000.0)
	while not end:
		changes, only_abs = em_update(args, total_bases, uniq_prop,
			only_abs, taxid2info, multimapped)
		if args.verbose:
			em_iter += 1
			echo('EM iteration ' + str(em_iter
				) + ' -- Sum of changes in abundances:' + str(
				sum([abs(v) for k,v in changes.items()])))
		end, prev_magnitude = end_condition(changes, prev_magnitude)

	for taxid in only_abs:  # incorporate updated abundance estimates
		taxids2abs[taxid][1] = only_abs[taxid]
	return taxids2abs


# Assign multimapped reads to specific organisms via proportional method,
#  	e.g. proportional to uniquely mapped reads; method used by MiCoP1
def resolve_multi_prop(args, taxids2abs, multimapped, taxid2info):
	# ensures early read assignments don't affect proportions of later ones
	to_add = {}
	for read_hits in multimapped:
		# get abundances of all taxids in read_hits
		all_taxids = list(set([hit[0] for hit in read_hits
			if hit[0] in taxids2abs]))
		if len(all_taxids) == 0:  # all hits were to taxids with no unique hits
			continue
		taxid_abs = [taxids2abs[tax][1] for tax in all_taxids]

		## now get cumulative proportional abs. of taxids relative to each other
		sumabs = sum(taxid_abs)
		if sumabs == 0.0:
			continue
		proportions = [ab / sumabs for ab in taxid_abs]
		hitlen = read_hits[0][-1]

		# divide hit length proportionally among hit taxids; divided assignment
		for i in range(len(all_taxids)):
			this_hitlen = proportions[i] * hitlen
			if not args.no_len_normalization:
				this_hitlen /= taxid2info[all_taxids[i]][0]
			if all_taxids[i] in to_add:
				to_add[all_taxids[i]] += this_hitlen
			else:
				to_add[all_taxids[i]] = this_hitlen


		#cumulative = [sum(proportions[:i+1]) for i in range(len(proportions))]

		# randomly proportionally choose which taxid to assign hit to
		#rand_draw = random.random()
		#for i in range(len(cumulative)):
		#	if rand_draw < cumulative[i] or i+1 == len(cumulative):
		#		assigned_taxid = all_taxids[i]
		#assigned_hits = [hit for hit in read_hits if hit[0] == assigned_taxid]

		# set the amount to add to taxid abundances
		#hitlen = assigned_hits[0][1]
		#if len(assigned_hits) == 2:  # paired end
		#	hitlen += assigned_hits[1][1]
		#if not args.no_len_normalization:
		#	hitlen /= taxid2info[assigned_taxid][0]
		#if assigned_taxid in to_add:
		#	to_add[assigned_taxid] += hitlen
		#else:
		#	to_add[assigned_taxid] = hitlen

	for taxid in to_add:  # add in the multimapped portions
		taxlin = taxid2info[taxid][-1].split('|')
		for tax in taxlin:
			if tax != '':
				taxids2abs[tax][1] += to_add[taxid]
		#taxids2abs[taxid][1] += to_add[taxid]
	return taxids2abs


# Renormalize each taxonomic rank so each rank sums to 100% abundance
def rank_renormalize(args, clades2abs, only_strains=False):
	rank_totals = {i:0.0 for i in RANKS}  # current rank abundance sums
	mapped_pct = 100.0
	if args.quantify_unmapped:  # only normalize against the pct of mapped reads
		mapped_pct = 100.0 - (100.0 * clades2abs['Unmapped'][-1])
	for clade in clades2abs:
		if clade == 'Unmapped':
			continue
		rank, ab = clades2abs[clade][1], clades2abs[clade][-1]
		if only_strains and rank != 'strain':
			continue
		rank_totals[rank] += ab  # add this to the rank sum total

	for clade in clades2abs:  # here's the normalization
		if clade == 'Unmapped':
			continue
		rank = clades2abs[clade][1]
		if only_strains and rank != 'strain':
			continue
		clades2abs[clade][-1]/= (rank_totals[clades2abs[clade][1]] / mapped_pct)
	return clades2abs


# given initital taxids2abs, some taxids will be higher than strain level;
#  	we insert "unknown" taxa down to strain level for normalization purposes
def gen_lower_taxa(taxids2abs):
	to_add = {}  # lower taxa to add after iteration
	for taxid in taxids2abs:
		taxid, rank, taxlin, namelin, ab = taxids2abs[taxid]
		if rank == 'strain':  # already lowest level
			continue
		# make a new taxid and name for this strain indicating that it is a
		#  	placeholder for a higher taxon that we cannot further narrow down
		rankpos = RANKS.index(rank)
		lowest_name = namelin.split('|')[rankpos]
		new_name = lowest_name + ' unknown strain'
		new_taxid = taxid+'.0'
		to_add[new_taxid] = [new_taxid, 'strain', taxlin + new_taxid,
			namelin + new_name, ab]

	# add in the new taxa
	for taxa in to_add:
		taxids2abs[taxa] = to_add[taxa]
	# clean out higher taxa listings -- they will return when we fill the tree
	#taxids2abs = {k:v for k,v in taxids2abs.items() if v[1] == 'strain'}
	return taxids2abs


# Get abundances for all clades in the tree and put in CAMI format
def tree_results_cami(args, taxids2abs):
	# rearrange fields to be in CAMI format
	for taxid in taxids2abs:
		old = taxids2abs[taxid]
		taxids2abs[taxid] = [taxid, old[3], old[5], old[4], old[1]]
	taxids2abs = gen_lower_taxa(taxids2abs)  # ensures everything strain level
	# always renormalize strains, to ensure legitimate profile
	#taxids2abs = rank_renormalize(args, taxids2abs, only_strains=True)

	# Now compute higher clade abundances
	clades2abs = {k:v for k,v in taxids2abs.items()}
	'''for taxid in taxids2abs:
		taxlin = taxids2abs[taxid][2].split('|')
		namelin = taxids2abs[taxid][3].split('|')
		for i in range(len(taxlin)-1):
			clade = taxlin[i]
			if clade == '':  # unspecified at this level
				continue
			if clade in clades2abs:  # already have clade entry, add abundance
				clades2abs[clade][-1] += taxids2abs[taxid][-1]
			else:
				# determine CAMI fields for this clade
				clade_taxid, clade_rank = clade, RANKS[i]
				clade_taxlin = '|'.join(taxlin[:i+1])
				clade_namelin = '|'.join(namelin[:i+1])
				clade_ab = taxids2abs[taxid][-1]  # currently just lower tax ab
				clades2abs[clade_taxid] = [clade_taxid, clade_rank,
					clade_taxlin, clade_namelin, clade_ab]'''

	if not args.no_rank_renormalization:
		clades2abs = rank_renormalize(args, clades2abs)
	return clades2abs


# Apply read cutoff and min_abundance thresholds, but latter relative to parent taxa
def prune_tree(args, taxids2abs):
	del_list = []
	for taxid in taxids2abs:
		if taxid == 'Unmapped':
			continue
		taxinfo = taxids2abs[taxid]
		num_reads, taxlin = taxinfo[0], taxinfo[-1]
		num_reads = taxinfo[0]
		if num_reads < args.read_cutoff:
			del_list.append(taxid)
			continue
		taxlin_splits = taxlin.split('|')
		if len(taxlin_splits) == 1:  # no parent
			continue
		for i in range(2, len(taxlin_splits) + 1):
			parent = taxlin.split('|')[-i]
			if parent != '':
				break
		parent_reads = taxids2abs[parent][0]
		if float(num_reads) / float(parent_reads) < args.min_abundance:
			del_list.append(taxid)

	# ensure the children of pruned parents are also pruned
	for taxid in taxids2abs:
		taxlin = taxids2abs[taxid][-1].split('|')
		for tax in taxlin:
			if tax in del_list:
				del_list.append(taxid)
	del_list = list(set(del_list))
	for taxid in del_list:
		del taxids2abs[taxid]
	return taxids2abs


# Processes uniquely-mapped reads, then estimates abundances using
#  	uniquely-mapped abundances and multimapped read information
def compute_abundances(args, infile, acc2info, tax2info):
	# taxids and higher clades to abundances, and multimapped reads dict
	taxids2abs, clades2abs, multimapped = {}, {}, {}
	if args.verbose:
		echo('Reading input file ' + infile)
	# run mapping and process to get uniq map abundances & multimapped reads
	taxids2abs, multimapped, total_bases, mmap_bases = map_and_process(
		args, infile, acc2info, tax2info)
	if args.verbose:
		echo('Done reading input file ' + infile)
	# filter out organisms below the read cutoff set by the user, if applicable
	#taxids2abs = {k:v for k,v in taxids2abs.items() if v[0] > args.read_cutoff}
	taxids2abs = prune_tree(args, taxids2abs)

	if args.verbose and args.assignment != 'none':
		echo('Assigning multimapped reads...')
	if len(multimapped) > 0 and args.assignment == 'em':
		taxids2abs = resolve_multi_em(args, total_bases, mmap_bases, taxids2abs,
			multimapped, tax2info)
	elif len(multimapped) > 0 and args.assignment == 'proportional':
		taxids2abs = resolve_multi_prop(args, taxids2abs, multimapped, tax2info)
	if args.verbose and args.assignment != 'none':
		echo('Multimapped reads assigned.')

	results = tree_results_cami(args, taxids2abs)
	if args.verbose:
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
		for clade in file_res:
			if clade not in results:
				results[clade] = file_res[clade]
			else:
				results[clade][-1] += file_res[clade][-1]

	if args.verbose:
		echo('Compiling and writing results...')
	rank_results = {i:[] for i in range(len(RANKS))}
	for clade in results:
		ab = results[clade][-1]  # avg over input files, truncate to 4 digits
		results[clade][-1] = float(math.trunc(
			(ab / len(args.infiles)) * (10**4))) / (10**4)
		rank = RANKS.index(results[clade][1])
		if rank == 7:  # strain; add extra CAMI genomeID and OTU fields
			taxid = results[clade][0]
			cami_genid, cami_otu = taxid, taxid.split('.')[0]
			results[clade].extend([cami_genid, cami_otu])
		rank_results[rank].append(results[clade])  # results per rank
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
			if not args.strain_level and i == len(RANKS)-1:
				continue  # skip the strain level unless user wants it
			lines = rank_results[i]  # all lines to write for this tax level
			# now sort clades in rank by descending abundance
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
	if args.pct_id == -1:  # not set by user
		args.pct_id = 0.5
		#if args.assignment == 'em':  # em assignment default
		#	args.pct_id = 0.5
		#else:  # proportional default
		#	args.pct_id = 0.95
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
