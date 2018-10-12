import argparse, operator, os, random, sys, time


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
		default='proportional', help='Method for assignming multimapped reads.')
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
	acc2info, taxid2info = {}, {}
	with(open(db_info, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			acclen = int(acclen)
			acc2info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in taxid2info:
				taxid2info[taxid][0] += acclen
			else:
				rank = get_taxid_rank(taxlin)
				taxid2info[taxid] = [acclen, rank, namelin, taxlin]
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
	sndary = (flag & 256 != 0)  # secondary hit
	# "bad" here means unmapped, chimeric mapping, or cigar string unavailable
	is_bad = (flag & 4 != 0) or (flag & 2048 != 0) or (cigar == '*')
	return pair1, pair2, sndary, is_bad


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


# Given all hits for a read, decide whether to use it and whether multimapped
# Returns multimapped reads, and if not, taxid and hit length of uniq map
def process_read(read_hits, pair1, pair2, pair1maps, pair2maps):
	if pair1 or pair2:  # paired read
		if pair1maps + pair2maps == 1:
			# if uniq mapped to one end and unmapped to other, count mapped end
			return [], read_hits[0][2], int(len(read_hits[0][9]))

		intersect, intersect_hits = intersect_read_hits(
											read_hits, pair1maps, pair2maps)
		if len(intersect) == 0:  # one end unmapped, other multimapped
			return [], 'Ambiguous', -1  # we consider this case too ambiguous
		elif len(intersect) == 1:  # read pairs only agree on one taxid
			hitlen = int(len(read_hits[0][9])) + int(len(read_hits[1][9]))
			return [], read_hits[0][2], hitlen  # we consider this a uniq map
		else:
			return intersect_hits, '', 0  #both pairs multimapped, process later

	else:  # single end
		if pair1maps > 1:  # multimapped
			return read_hits, '', 0  # multimapped, acc and hitlen not needed
		else:
			return False, read_hits[0][2], int(len(read_hits[0][9]))


# Reads sam file -> returns a dict of taxids mapped to number of uniquely mapped
# 	reads and bases for that taxid, and a list of multimapped reads.
def process_samfile(args, samfile, acc2info, taxid2info):
	taxids2abs, multimapped = {}, []  # taxids to abundances, multimapped reads
	prev_read, read_hits = '', [] # read tracker, all full hits for read pair
	pair1maps, pair2maps = 0, 0  # reads mapped to each pair (single = pair1)

	with(open(samfile, 'r')) as infile:
		if args.verbose:
			print('Reading sam file ' + samfile)
		for line in infile:
			if line.startswith('@'):
				continue

			splits = line.strip().split()
			pair1, pair2, sndary, is_bad = parse_flag(int(splits[1]), splits[5])
			if is_bad:  #unmapped, chimeric mapping, or cigar string unavailable
				continue
			if filter_line(args, splits):
				continue  # filter read based on user specifications
			# here we change accession to taxid since we want to
			#  	profile by organisms, not accessions
			splits[2] = acc2info[splits[2]][1]

			read, ref = splits[0], splits[2]
			if read != prev_read:
				# get uniq hit taxid and hitlen, or multimapped hits intersect
				intersect_hits, taxid, hitlen = process_read(
					read_hits, pair1, pair2, pair1maps, pair2maps)
				# reset these read-specific variables
				prev_read, read_hits, pair1maps, pair2maps = read, [splits], 0,0
				if taxid == 'Ambiguous':  # ambiguous mapping (see process_read)
					continue
				if taxid != '' and not args.no_len_normalization:
					hitlen /= taxid2info[taxid][0]  # normalize by genome length
				if intersect_hits == []:  # unique hit
					if taxid in taxids2abs:
						taxids2abs[taxid][0] += 1  # reads hit
						taxids2abs[taxid][1] += hitlen  # bases hit
					else:  # also store lineage for taxid
						taxids2abs[taxid] = [1, hitlen] + taxid2info[taxid]
				else:  # multimapped hit
					multimapped.append(intersect_hits)
			else:
				pair1maps += pair1 or not(pair1 or pair2)  # pair1 or single
				pair2maps += pair2  # unchanged if pair2 false
				read_hits.append(splits)

		if args.verbose:
			print('Done reading sam file ' + samfile)
	return taxids2abs, multimapped


# Use EM algorithm to revise taxid abundances using multimapped information,
#  	without explicitly assigning reads. Recommended method in MiCoP2.
def resolve_multi_em(args, taxids2abs, multimapped):
	pass  # to be implemented


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
		if only_strains and clades2abs[clade] != 'strain':
			continue
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
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[2] != 'strain'}
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
		for i in range(len(taxlin)):
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
	taxids2abs, multimapped = process_samfile(args, samfile, acc2info, tax2info)
	if args.verbose:
		print('Done reading sam file ' + samfile)
	# filter out organisms below the read cutoff set by the user, if applicable
	taxids2abs = {k:v for k,v in taxids2abs.items() if v[0] > args.read_cutoff}

	if args.verbose:
		print('Assigning multimapped reads...')
	if args.assignment == 'em':
		taxids2abs = resolve_multi_em(args, taxids2abs, multimapped)
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
		results[clade][-1] /= len(args.sam)  # average over all input files
		rank = RANKS.index(results[clade][1])
		if rank == 7:  # strain; add extra CAMI genomeID and OTU fields
			taxid = results[clade][0]
			if taxid.endswith('.0'):  # unidentified strain
				cami_genid, cami_otu = taxid, taxid.split('.')[0]
			else:
				cami_genid, cami_otu = taxid + '.1', taxid
			results[clade].extend([cami_genid, cami_otu])
		rank_results[rank].append(results[clade])  # results per rank
	return rank_results


def write_results(args, rank_results):
	with(open(args.output, 'w')) as outfile:
		# Print some CAMI format header lines
		outfile.write('@SampleID:' + ','.join(args.sam) + '\n')
		outfile.write('@Version:MiCoP2-v0.1\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\
			\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

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
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')


def main():
	args = parseargs()
	if args.pct_id == -1:  # not set by user
		if args.assignment == 'em':  # em assignment default
			args.pct_id = 0.0
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
