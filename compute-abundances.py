import argparse, math, os, random, sys


__location__ = os.path.realpath(os.path.join(os.getcwd(),
								os.path.dirname(__file__)))
db_info = __location__ + '/data/db_info.txt'
RANKS = ['superkingdom', 'phylum', 'class', 'order',
		'family', 'genus', 'species', 'strain']


def parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Compute abundance estimations for species in a sample.')
	parser.add_argument('sam', nargs='+',
		help='.sam file or list of .sam files to process. Required.')
	parser.add_argument('--assignment', choices=['em', 'proportional'],
		default='em', help='Method for assignming multimapped reads.')
	parser.add_argument('--min_abundance', type=float, default=10**-10,
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
	parser.add_argument('--output', default='abundances.txt',
		help='Output abundances file. Default: abundances.txt')
	parser.add_argument('--pct_id', type=float, default=-1,
		help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--read_cutoff', type=int, default=-1,
		help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to sam file name(s).')
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


# Parses information in db_info.txt file, which maps NCBI accession to
# 	accession length, taxid, name lineage, and taxid lineage.
# Also maps all lowest-level taxids to organism length, which is the sum of
#  	accession lengths for accessions with that taxid, and lineage info
def get_acc2info(args):
	if args.verbose:
		print('Reading db_info.txt...')
	acc2info, taxid2info = {}, {}
	with(open(db_info, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			rank = get_taxid_rank(taxlin)
			if rank == 'strain':
				taxid += '.1'  # CAMI formatting specification
				taxlin += '.1'
			acclen = int(acclen)
			acc2info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in taxid2info:
				taxid2info[taxid][0] += acclen
			else:
				taxid2info[taxid] = [acclen, rank, namelin, taxlin]
	if args.verbose:
		print('Done reading db_info.txt.')
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


# Remove hits to taxids not in taxids2abs (no unique mappings to that taxid)
def preprocess_multimapped(total_bases, multimapped, taxids2abs):
	for i in range(len(multimapped)):
		readquals, paired = multimapped[i][0][10], multimapped[i][0][-1]
		multimapped[i] = [hit for hit in multimapped[i] if hit[2] in taxids2abs]
		if len(multimapped[i]) > 0:
			multimapped[i][0][10], multimapped[i][0][-1] = readquals, paired
			total_bases += len(readquals)
	multimapped = [read for read in multimapped if len(read) > 0]
	return total_bases, multimapped


# Reads sam file -> returns a dict of taxids mapped to number of uniquely mapped
# 	reads and bases for that taxid, and a list of multimapped reads.
def process_samfile(args, samfile, acc2info, taxid2info):
	taxids2abs, multimapped = {}, []  # taxids to abundances, multimapped reads
	prev_read, read_hits = '', [] # read tracker, all hits for read (full lines)
	pair1maps, pair2maps = 0, 0  # reads mapped to each pair (single = pair1)
	total_bases = 0

	with(open(samfile, 'r')) as infile:
		for line in infile:
			if line.startswith('@'):
				continue
			splits = line.strip().split()
			pair1, pair2, chimer, is_bad = parse_flag(int(splits[1]), splits[5])
			if is_bad:  # unmapped or cigar string unavailable
				continue
			splits[1] = [pair1, pair2, chimer, is_bad]  # store flag fields
			# here we change accession to taxid since we want to
			#  	profile by organisms, not accessions
			splits[2] = acc2info[splits[2]][1]

			read, ref = splits[0], splits[2]
			if read != prev_read:
				# get uniq hit taxid and hitlen, or multimapped hits intersect
				intersect_hits, taxid, readquals, hitlen = process_read(
					args, read_hits, pair1, pair2, pair1maps, pair2maps)
				# reset these read-specific variables
				prev_read, read_hits, pair1maps, pair2maps = read, [splits], 0,0
				if taxid == 'Ambiguous':  # ambiguous mapping (see process_read)
					continue
				if taxid != '' and not (args.no_len_normalization
											or args.assignment == 'em'):
					hitlen /= taxid2info[taxid][0]  # normalize by genome length
				if intersect_hits == []:  # unique hit
					total_bases += hitlen
					if taxid in taxids2abs:
						taxids2abs[taxid][0] += 1  # reads hit
						taxids2abs[taxid][1] += hitlen  # bases hit
					else:  # also store lineage for taxid
						taxids2abs[taxid] = [1, hitlen] + taxid2info[taxid]
				else:  # multimapped hit
					intersect_hits[0][10] = readquals  # ensure qual scores
					intersect_hits[0].append(pair1 or pair2)  # paired or not
					multimapped.append(intersect_hits)
			else:
				pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
				pair2maps += pair2  # unchanged if pair2 false
				read_hits.append(splits)
	total_bases, multimapped = preprocess_multimapped(total_bases, multimapped,
															taxids2abs)
	return taxids2abs, multimapped, float(total_bases)


# Initial EM abundances estimate: similar to Karp paper, except proportion of
#  	uniquely-mapped reads is changed to (possibly length-normalized) proportion
#  	of uniquely-mapped bases, allowing variable read length to factor in
def initital_estimate(taxids2abs, multimapped, total_bases):
	# calculate share of total bases uniquely mapped to each taxid
	taxids2abs = {k: v/total_bases for k,v in taxids2abs.items()}
	uniq_proportions = {k:v for k,v in taxids2abs.items()}

	# now add in initial estimate for multimapped reads: evenly-distributed
	all_multimap_len = float(sum([len(read[0][10]) for read in multimapped]))
	all_mmap_taxids = []  # set of multimapped taxids with >= 1 uniq mapped read
	for read in multimapped:
		all_mmap_taxids.extend([hit[2] for hit in read if hit[2] in taxids2abs])
	all_mmap_taxids = set(all_mmap_taxids)
	even_dist = float(all_multimap_len) / float(total_bases) / float(len(all_mmap_taxids))
	for taxid in taxids2abs:
		taxids2abs[taxid] += even_dist
	return uniq_proportions, taxids2abs


# Based on magnitude of changes since last EM estimate, decide whether to stop
def end_condition(changes):
	ab_magnitude = sum([abs(v) for k,v in changes.items()])
	if ab_magnitude < (10 ** -3):
		return True
	return False


def single_read_likelihood(probs, cigar, taxid_ab):
	# curval = current cigar letter amount, cur_ind = base probs. index
	base_likelihoods, curval, cur_ind = [], 0, 0
	for ch in cigar:
		if ch.isdigit():
			curval = (curval * 10) + int(ch)
		else:
			if ch == 'M' or ch == '=':  # base likelihoods for matched bases
				base_likelihoods.extend(
					[1-probs[i] for i in range(cur_ind, cur_ind + curval)])
			else:  # base likelihoods for mismatched base
				base_likelihoods.extend(
					[probs[i]/3.0 for i in range(cur_ind, cur_ind + curval)])
			cur_ind = cur_ind + curval
			curval = 0
	# take product of base_likelihoods, times previous EM taxid abundance est.
	result = taxid_ab
	for bl in base_likelihoods:
		result *= bl
	return result


# computes likelihood of read being assigned to one hit, divided over the total
#  	likelihood of all hit assignments
def read_likelihood_ratio(read_hits, taxids2abs):
	base_quals = read_hits[0][10]
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
		likelihoods[taxid] = single_read_likelihood(prob_bases_wrong, cigars,
														taxids2abs[taxid])

	# compute likelihoods as proportions of sum over all assignment likelihoods
	total_likelihood = sum([v for k, v in likelihoods.items()])
	likelihoods = {k: v / total_likelihood for k,v in likelihoods.items()}
	return likelihoods


def em_update(args, total_bases, uniq_prop, taxids2abs,taxid2info, multimapped):
	init_abs = {k:v for k,v in taxids2abs.items()}
	# get read assignment likelihoods for each read
	likelihood_sums = {}
	for read in multimapped:
		# base_share is the share of total bases reflected in this read
		base_share = float(len(read[0][10])) / float(total_bases)
		if read[0][-1]:  # if paired end, account for other end
			base_share += float(len(read[-1][9])) / float(total_bases)
		read_likelihoods = read_likelihood_ratio(read, taxids2abs)
		for taxid in read_likelihoods:
			if taxid not in likelihood_sums:
				likelihood_sums[taxid] = (read_likelihoods[taxid] * base_share)
			else:
				likelihood_sums[taxid] += (read_likelihoods[taxid] * base_share)

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


def min_likelihood_filter(multimapped, taxids2abs):
	for read in range(len(multimapped)):
		del_hits = []  # hits to be deleted
		base_quals, paired = multimapped[read][0][10], multimapped[read][0][-1]
		hit_dict = {}  # group reads by reference hit, needed for paired end
		for hit in multimapped[read]:
			if hit[2] in hit_dict:
				hit_dict[hit[2]].append(hit)
			else:
				hit_dict[hit[2]] = [hit]

		# compute read likelihoods
		for taxid in hit_dict:
			cigars = hit_dict[taxid][0][5]
			if len(hit_dict[taxid]) == 2:  # paired end
				cigars += hit_dict[taxid][1][5]
			else:
				cigars += hit_dict[taxid][0][5]  # see: read_likelihood_ratio
			# compute prob. of each base being wrong, then this read's likelihood
			prob_bases_wrong = [10 ** (-(ord(ch) - 33)/10) for ch in base_quals]
			likelihood = single_read_likelihood(prob_bases_wrong, cigars,
															taxids2abs[taxid])
			if likelihood  < 10 ** -300:  # min. likelihood
				del_hits.extend([hit for hit in hit_dict[taxid]])

		# clean extremely low-likelihood hits and reads with only these bad hits
		multimapped[read] = [i for i in multimapped[read] if i not in del_hits]
		if len(multimapped[read]) > 0:
			multimapped[read][0][10],multimapped[read][0][-1]=base_quals,paired
	multimapped = [read for read in multimapped if len(read) > 0]
	return multimapped


# Use EM algorithm to revise taxid abundances using multimapped information,
#  	without explicitly assigning reads. Recommended method in MiCoP2.
def resolve_multi_em(args, total_bases, taxids2abs, multimapped, taxid2info):
	# we temporarily just keep abundance to simplify EM operations
	only_abs = {k: v[1] for k,v in taxids2abs.items()}
	# filter out hits and reads with extremely low likelihood
	multimapped = min_likelihood_filter(multimapped, only_abs)
	if len(multimapped) == 0:
		return taxids2abs
	uniq_prop, only_abs = initital_estimate(only_abs, multimapped, total_bases)

	# we have a dict showing changes in ab. estimates for each step;
	#  	 if this changes dict has small enough magnitude, end EM updates
	changes = {k:(only_abs[k]-v) for k,v in uniq_prop.items()}
	iter = 0
	while not end_condition(changes):
		changes, only_abs = em_update(args, total_bases, uniq_prop,
										only_abs, taxid2info, multimapped)
		if args.verbose:
			iter += 1
			print('EM iteration', str(iter), '-- Sum of changes in abundances:',
			 	str(sum([abs(v) for k,v in changes.items()])))

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
		all_taxids = list(set(
						[hit[2] for hit in read_hits if hit[2] in taxids2abs]))
		if len(all_taxids) == 0:  # all hits were to taxids with no unique hits
			continue
		taxid_abs = [taxids2abs[tax][1] for tax in all_taxids]

		# now get cumulative proportional abs. of taxids relative to each other
		sumabs = sum(taxid_abs)
		proportions = [ab / sumabs for ab in taxid_abs]
		cumulative = [sum(proportions[:i+1]) for i in range(len(proportions))]

		# randomly proportionally choose which taxid to assign hit to
		rand_draw = random.random()
		for i in range(len(cumulative)):
			if rand_draw < cumulative[i] or i+1 == len(cumulative):
				assigned_taxid = all_taxids[i]
		assigned_hits = [hit for hit in read_hits if hit[2] == assigned_taxid]

		# set the amount to add to taxid abundances
		hitlen = len(assigned_hits[0][9]) + len(assigned_hits[1][9])
		if not args.no_len_normalization:
			hitlen /= taxid2info[assigned_taxid][0]
		if assigned_taxid in to_add:
			to_add[assigned_taxid] += hitlen
		else:
			to_add[assigned_taxid] = hitlen

	for taxid in to_add:  # add in the multimapped portions
		taxids2abs[taxid][1] += to_add[taxid]
	return taxids2abs


# Renormalize each taxonomic rank so each rank sums to 100% abundance
def rank_renormalize(clades2abs, only_strains=False):
	rank_totals = {i:0.0 for i in RANKS}  # current rank abundance sums
	for clade in clades2abs:
		rank, ab = clades2abs[clade][1], clades2abs[clade][-1]
		rank_totals[rank] += ab  # add this to the rank sum total
	for clade in clades2abs:  # here's the normalization
		clades2abs[clade][-1] /= (rank_totals[clades2abs[clade][1]] / 100.0)
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
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[1] == 'strain'}
	return taxids2abs


# Get abundances for all clades in the tree and put in CAMI format
def tree_results_cami(args, taxids2abs):
	# rearrange fields to be in CAMI format
	for taxid in taxids2abs:
		old = taxids2abs[taxid]
		taxids2abs[taxid] = [taxid, old[3], old[5], old[4], old[1]]
	taxids2abs = gen_lower_taxa(taxids2abs)  # ensures everything strain level
	# always renormalize strains, to ensure legitimate profile
	taxids2abs = rank_renormalize(taxids2abs, only_strains=True)

	# Now compute higher clade abundances
	clades2abs = {k:v for k,v in taxids2abs.items()}
	for taxid in taxids2abs:
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
										  clade_taxlin, clade_namelin, clade_ab]

	if not args.no_rank_renormalization:
		clades2abs = rank_renormalize(clades2abs)
	return clades2abs


def compute_abundances(args, samfile, acc2info, tax2info):
	# taxids and higher clades to abundances, and multimapped reads dict
	taxids2abs, clades2abs, multimapped = {}, {}, {}
	if args.verbose:
		print('Reading sam file ' + samfile)
	# parse sam file to get uniq map abundances and multimapped reads
	taxids2abs, multimapped, total_bases = process_samfile(args, samfile,
															acc2info, tax2info)
	if args.verbose:
		print('Done reading sam file ' + samfile)
	# filter out organisms below the read cutoff set by the user, if applicable
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[0] > args.read_cutoff}

	if args.verbose:
		print('Assigning multimapped reads...')
	if args.assignment == 'em':
		taxids2abs = resolve_multi_em(args, total_bases, taxids2abs,
										multimapped, tax2info)
	else:  # proportional
		taxids2abs = resolve_multi_prop(args, taxids2abs, multimapped, tax2info)
	if args.verbose:
		print('Multimapped reads assigned.')

	results = tree_results_cami(args, taxids2abs)
	if args.verbose:
		print('Done computing abundances for sam file ' + samfile)
	return results


def gather_results(args, acc2info, taxid2info):
	results = {}
	for sam in args.sam:
		if args.verbose:
			print('Computing abundances for sam file ' + sam + '...')
		sam_res = compute_abundances(args, sam, acc2info, taxid2info)
		for clade in sam_res:
			if clade not in results:
				results[clade] = sam_res[clade]
			else:
				results[clade][-1] += sam_res[clade][-1]

	if args.verbose:
		print('Compiling and writing results...')
	rank_results = {i:[] for i in range(len(RANKS))}
	for clade in results:
		ab = results[clade][-1]  # avg over all sam files, truncate to 4 digits
		results[clade][-1] = float(math.trunc(
								(ab / len(args.sam)) * (10 ** 4))) / (10 ** 4)
		rank = RANKS.index(results[clade][1])
		if rank == 7:  # strain; add extra CAMI genomeID and OTU fields
			taxid = results[clade][0]
			cami_genid, cami_otu = taxid, taxid.split('.')[0]
			results[clade].extend([cami_genid, cami_otu])
		rank_results[rank].append(results[clade])  # results per rank
	return rank_results


def write_results(args, rank_results):
	with(open(args.output, 'w')) as outfile:
		# Print some CAMI format header lines
		if args.sampleID == 'NONE':
			outfile.write('@SampleID:' + ','.join(args.sam) + '\n')
		else:
			outfile.write('@SampleID:' + args.sampleID + '\n')
		outfile.write('@Version:MiCoP2-v0.1\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t' +
			'PERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

		for i in range(len(RANKS)):
			lines = rank_results[i]  # all lines to write for this tax level
			# now sort clades in rank by descending abundance
			if i != 7:  # not strains, which have extra fields
				lines.sort(key=lambda x: 100.0-x[-1])
			else:  # strains
				lines.sort(key=lambda x: 100.0-x[-3])
			if lines == None or len(lines) < 1:
				continue
			for line in lines:
				if line[4] < args.min_abundance:
					continue
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')


def main():
	args = parseargs()
	if args.pct_id == -1:  # not set by user
		if args.assignment == 'em':  # em assignment default
			args.pct_id = 0.5
		else:  # proportional default
			args.pct_id = 0.95
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		print('Error: --pct_id must be between 0.0 and 1.0, inclusive.')
		sys.exit()
	open(args.output, 'w').close()  # test to see if writeable

	# maps NCBI accession to length, taxid, name lineage, taxid lineage
	acc2info, taxid2info = get_acc2info(args)
	# gathers results for all sam files, combines, and organizes into tax levels
	rank_results = gather_results(args, acc2info, taxid2info)
	write_results(args, rank_results)


if __name__ == '__main__':
	main()
#
