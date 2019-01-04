import argparse, math, os, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Run CMash and" +
				" select a subset of the whole database to align to.")
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('db_dir',
		help='Directory with all organism files in the full database.')
	parser.add_argument('--cutoff', type=float, default=-1.0,
		help = 'CMash cutoff value. Default is 1/(log10(reads file bytes)**2).')
	parser.add_argument('--output', default='cmashed_db.fna',
		help = 'Name of directory to write selected genome files to.')
	args = parser.parse_args()
	return args


def run_cmash_and_cutoff(args):
	cmash_db_loc = 'data/cmash_db.h5'
	# temporarily write cmash results to the --output file; should be writeable
	cmash_proc = subprocess.Popen(['CMash/scripts/QueryDNADatabase.py',
		args.reads, cmash_db_loc, args.output]).wait()

	organisms_to_include = []
	with(open(args.output, 'r')) as cmash_results:
		cmash_results.readline()  # skip header line
		for line in cmash_results:
			splits = line.strip().split(',')
			organism, containment_index = splits[0], float(splits[-1])
			if containment_index >= args.cutoff:
				organisms_to_include.append(organism)
	return organisms_to_include


def make_db(args, organisms_to_include):
	open(args.output, 'w').close()  # clear cmash results; no longer needed
	with(open(args.output, 'a')) as outfile:
		for organism in organisms_to_include:
			organism_fname = args.db_dir + organism
			# write organisms to full db via cat to append-mode file handler
			subprocess.call(['cat', organism_fname], stdout=outfile)


def main():
	args = parseargs()
	if args.cutoff == -1.0:  # not set by user
		fsize_bytes = os.path.getsize(args.reads)
		args.cutoff = 1.0 / (math.log10(fsize_bytes) ** 2)
	elif args.cutoff < 0.0 or args.cutoff > 1.0:
		print('Error: args.cutoff must be between 0 and 1, inclusive.')
		sys.exit()
	if not args.db_dir.endswith('/'):
		args.db_dir += '/'

	organisms_to_include = run_cmash_and_cutoff(args)
	make_db(args, organisms_to_include)


if __name__ == '__main__':
	main()
#
