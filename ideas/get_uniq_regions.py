import argparse, glob, multiprocessing, random, subprocess, sys, time


start = time.time()  # start a program timer
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
	parser = argparse.ArgumentParser(description='Get unique regions of accessions.')
	parser.add_argument('db_info', help='db_info file to get taxonomic info from.')
	parser.add_argument('organism_files',
		help='Directory with individual organism genome files.')
	parser.add_argument('temp_dir',
		help='Temporary working directory (must already exist).')
	parser.add_argument('--max_neighbors', type=int, default=10,
		help='Max neighbors for each taxid to run nucmer against.')
	parser.add_argument('--nucmer_files_exist', action='store_true',
		help='Use this flag if nucmer output files already exist.')
	parser.add_argument('--num_procs', type=int, default=10,
		help='Size of multiprocessing pool.')
	parser.add_argument('--output', default='uniq_regions.txt',
		help='Output file. Default: uniq_regions.txt')
	parser.add_argument('--print_only', action='store_true',
		help='Instead of running MUMmer, just write the lines that would be run.')
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
def get_acc_to_info(dbinfo):
	acc_to_info, tax_to_info = {}, {}
	with(open(dbinfo, 'r')) as infofile:
		infofile.readline()  # skip header line
		for line in infofile:
			acc, acclen, taxid, namelin, taxlin = line.strip().split('\t')
			rank = get_taxid_rank(taxlin)
			#if rank == 'strain' and acc != 'Unmapped':
			#	taxid += '.1'  # CAMI formatting specification
			#	taxlin += '.1'
			acclen = int(acclen)
			acc_to_info[acc] = [acclen, taxid, namelin, taxlin]
			if taxid in tax_to_info:
				tax_to_info[taxid][0] += acclen
			else:
				tax_to_info[taxid] = [acclen, rank, namelin, taxlin]
	return acc_to_info, tax_to_info


# Maps higher level taxids to all lowest-level (leaf) taxids under it
def get_taxid_to_leaves(tax_to_info):
	tax_to_leaves = {}
	for leaf in tax_to_info:
		taxlin = tax_to_info[leaf][-1]
		for taxid in taxlin:
			if taxid not in tax_to_leaves:
				tax_to_leaves[taxid] = []
			tax_to_leaves[taxid].append(leaf)
	for taxid in tax_to_leaves:
		tax_to_leaves[taxid] = list(set(tax_to_leaves[taxid]))
	return tax_to_leaves


# Get nearest neighbors of the taxid up to max_neighbors
def get_neighbors(tax_to_info, tax_to_leaves, max_neighbors, nucmer_files_exist, temp_dir):
	tax_to_neighbors = {}
	if nucmer_files_exist:  # just generate from nucmer file names
		for fname in glob.glob(temp_dir + 'coords*'):
			taxid, neighbor = fname.split('/')[-1].split('-')[1:3]
			neighbor = neighbor.split('.txt')[0]
			if taxid not in tax_to_neighbors:
				tax_to_neighbors[taxid] = []
			tax_to_neighbors[taxid].append(neighbor)
		return tax_to_neighbors

	for leaf in tax_to_info:
		tax_to_neighbors[leaf] = []
		cur_neighbors = 0
		taxlin = tax_to_info[leaf][-1].split('|')
		# iterate backwards through list, adding neighbors up to max_neighbors
		for i in range(2, len(taxlin) + 1):
			taxid = taxlin[-i]
			if taxid == '' or taxid not in tax_to_leaves:
				continue
			leaf_neighbors = tax_to_leaves[taxid]
			leaf_neighbors.remove(leaf)  # don't measure leaf against itself
			if cur_neighbors + len(leaf_neighbors) <= max_neighbors:
				tax_to_neighbors[leaf].extend(leaf_neighbors)
				tax_to_neighbors[leaf] = list(set(tax_to_neighbors[leaf]))
				cur_neighbors = len(tax_to_neighbors[leaf])
			else:  # randomly select enough neighbors to hit max_neighbors
				random.shuffle(leaf_neighbors)
				amnt_to_add = max_neighbors - cur_neighbors
				leaf_neighbors = leaf_neighbors[:amnt_to_add]
				tax_to_neighbors[leaf].extend(leaf_neighbors)
				tax_to_neighbors[leaf] = list(set(tax_to_neighbors[leaf]))
				break
	return tax_to_neighbors


# Run nucmer and show-coords from MUMmer on a taxid-neighbor pair
def nucmer_show_coords_proc(org_files, temp_dir, taxid, neighbor):
	ref_file = org_files + 'taxid_' + taxid.replace('.', '_') + '_genomic.fna'
	query_file = org_files + 'taxid_' + neighbor.replace('.', '_') + '_genomic.fna'
	nucmer_outfile = temp_dir + 'nucmer-' + taxid + '-' + neighbor + '.txt'
	coords_outfile = temp_dir + 'coords-' + taxid + '-' + neighbor + '.txt'
	subprocess.Popen(['nucmer', ref_file, query_file, '--delta=' + nucmer_outfile]).wait()
	with(open(coords_outfile, 'w')) as proc_outfile:
		subprocess.Popen(['show-coords', nucmer_outfile], stdout = proc_outfile).wait()
	return 0


# Execute nucmer/show-coords on all taxid-neighbor pairs; parallelize via multiprocessing
def run_mummer_mp(tax_to_neighbors, org_files, temp_dir, num_procs):
	pool = multiprocessing.Pool(num_procs)
	for taxid in tax_to_neighbors:
		for neigh in tax_to_neighbors[taxid]:
			pool.apply_async(nucmer_show_coords_proc, (org_files, temp_dir, taxid, neigh,))
	pool.close()
	pool.join()


def print_runlines(tax_to_neighbors, org_files, temp_dir, output):
	with(open(args.output, 'w')) as outfile:
		for taxid in tax_to_neighbors:
			for neighbor in tax_to_neighbors[taxid]:
				ref_file = org_files + 'taxid_' + taxid.replace('.', '_') + '_genomic.fna'
				query_file = org_files + 'taxid_' + neighbor.replace('.', '_') + '_genomic.fna'
				nucmer_outfile = temp_dir + 'nucmer-' + taxid + '-' + neighbor + '.txt'
				coords_outfile = temp_dir + 'coords-' + taxid + '-' + neighbor + '.txt'
				nucmer_runline = ' '.join(['nucmer', ref_file, query_file, '--delta=' + nucmer_outfile])
				coords_runline = ' '.join(['show-coords', nucmer_outfile, '>', coords_outfile])
				outfile.write(nucmer_runline + '\n' + coords_runline + '\n')


# Parses show-coords output for each taxid-neighbor pair for a single taxid,
#    extracting the shared blocks for each accession
def get_acc_shared_blocks(taxid, temp_dir, tax_to_neighbors):
	acc_to_shared_blocks = {}
	#for fname in glob.glob(temp_dir + 'coords-' + taxid + '*'):
	for neighbor in tax_to_neighbors[taxid]:
		fname = temp_dir + 'coords-' + taxid + '-' + neighbor + '.txt'
		with(open(fname, 'r')) as coords_outfile:
			for i in range(5):
				coords_outfile.readline()  # skip headers
			for line in coords_outfile:
				splits = line.split()
				start, end, acc = int(splits[0]), int(splits[1]), splits[11]
				if 'gi|' in acc:
					acc = acc.split('|')[3]
				if acc not in acc_to_shared_blocks:
					acc_to_shared_blocks[acc] = []
				acc_to_shared_blocks[acc].append([start, end])
	return acc_to_shared_blocks


# Fuse shared blocks that overlap
def fuse_shared_blocks(acc_to_shared_blocks):
	fused_blocks = {}
	for acc in acc_to_shared_blocks:
		fused_blocks[acc] = []
		acc_to_shared_blocks[acc].sort(key = lambda x: x[0])
		for block in acc_to_shared_blocks[acc]:
			if len(fused_blocks[acc]) == 0:
				fused_blocks[acc].append(block)
				continue
			last_block_end = fused_blocks[acc][-1][1]
			# three cases: need new block, extend block, or don't need to extend
			if block[0] > last_block_end + 1:
				fused_blocks[acc].append(block)
			elif block[1] > last_block_end:
				fused_blocks[acc][-1][1] = block[1]
	return fused_blocks


# Return the unique blocks: the regions outside of the fused shared blocks
def get_acc_uniq_blocks(fused_blocks, acc_to_info):
	min_block_size = 100  # don't use unique blocks smaller than this
	acc_to_uniq_blocks = {}
	for acc in fused_blocks:
		acc_to_uniq_blocks[acc] = []
		acclen = acc_to_info[acc][0]
		if len(fused_blocks[acc]) == 0:  # no shared blocks; all unique
			acc_to_uniq_blocks[acc] = [[0, acclen]]
			continue
		if fused_blocks[acc][0][0] > 1:
			uniq_block_size = fused_blocks[acc][0][0]
			if uniq_block_size >= min_block_size:
				acc_to_uniq_blocks[acc].append([1, fused_blocks[acc][0][0]])

		# areas between shared blocks are unique
		for i in range(1, len(fused_blocks[acc])):
			prev_block_end = fused_blocks[acc][i - 1][1]
			next_block_start = fused_blocks[acc][i][0]
			uniq_block_size = next_block_start - prev_block_end
			if uniq_block_size < min_block_size:
				continue
			acc_to_uniq_blocks[acc].append([prev_block_end, next_block_start])

		if fused_blocks[acc][-1][1] < acclen:
			uniq_block_size = acclen - fused_blocks[acc][-1][1]
			if uniq_block_size >= min_block_size:
				acc_to_uniq_blocks[acc].append([fused_blocks[acc][-1][1], acclen])
	return acc_to_uniq_blocks


# For each taxid, compute the unique blocks for each of its accessions using the
#    show-coords output, then write to a file
def write_acc_uniq_blocks(output, tax_to_neighbors, acc_to_info, temp_dir):
	with(open(output, 'w')) as outfile:
		accs_written = {}
		for taxid in tax_to_neighbors:
			if taxid == 'Unmapped':
				continue
			acc_to_shared_blocks = get_acc_shared_blocks(taxid, temp_dir, tax_to_neighbors)
			fused_blocks = fuse_shared_blocks(acc_to_shared_blocks)
			acc_to_uniq_blocks = get_acc_uniq_blocks(fused_blocks, acc_to_info)
			for acc in acc_to_uniq_blocks:
				accs_written[acc] = True
				outfile.write(acc)
				for block in acc_to_uniq_blocks[acc]:
					outfile.write('\t' + str(block))
				outfile.write('\n')
		# write accs that are entirely unique and thus were not written above
		for acc in acc_to_info:
			if acc in accs_written:
				continue
			acclen = acc_to_info[acc][0]
			uniq_region = [1, acclen]  # whole genome is unique
			outfile.write(acc + '\t' + str(uniq_region) + '\n')


def main(args = None):
	if args == None:
		args = parseargs()
	if not args.organism_files.endswith('/'):
		args.organism_files += '/'
	if not args.temp_dir.endswith('/'):
		args.temp_dir += '/'

	echo('Reading db_info file...', args.verbose)
	acc_to_info, tax_to_info = get_acc_to_info(args.db_info)
	echo('Getting leaf nodes for each taxid...', args.verbose)
	tax_to_leaves = get_taxid_to_leaves(tax_to_info)
	echo('Getting neighbors for each leaf node...', args.verbose)
	tax_to_neighbors = get_neighbors(tax_to_info, tax_to_leaves,
		args.max_neighbors, args.nucmer_files_exist, args.temp_dir)
	if not args.print_only:
		if not args.nucmer_files_exist:
			echo('Running MUMmer processes...', args.verbose)
			run_mummer_mp(tax_to_neighbors, args.organism_files, args.temp_dir, args.num_procs)
		echo('Calculating unique blocks and writing output...', args.verbose)
		write_acc_uniq_blocks(args.output, tax_to_neighbors, acc_to_info, args.temp_dir)
	else:
		echo('Printing nucmer/show-coords run lines to --output...', args.verbose)
		print_runlines(tax_to_neighbors, args.organism_files, args.temp_dir, args.output)


if __name__ == '__main__':
	args = parseargs()
	main(args)
#
