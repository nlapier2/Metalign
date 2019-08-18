import argparse, math, subprocess, sys, time
# Import metalign modules
import select_db as select
import map_and_profile as mapper


def metalign_parseargs():  # handle user arguments
	parser = argparse.ArgumentParser(
		description='Runs full metalign pipeline on input reads file(s).')
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('data',
		help='Path to data/ directory with the files from setup_data.sh')
	parser.add_argument('--cutoff', type=float, default=-1.0,
		help='CMash cutoff value. Default is 1/(log10(reads file bytes)**2).')
	parser.add_argument('--db_dir', default = 'AUTO',
		help='Directory with all organism files in the full database.')
	parser.add_argument('--dbinfo_in', default='AUTO',
		help = 'Location of db_info file. Default: data/db_info.txt')
	parser.add_argument('--keep_temp_files', action = 'store_true',
		help='Keep KMC files instead of deleting after this script finishes.')
	parser.add_argument('--min_abundance', type=float, default=10**-4,
		help='Minimum abundance for a taxa to be included in the results.')
	parser.add_argument('--length_normalize', action='store_true',
		help='Normalize abundances by genome length.')
	parser.add_argument('--rank_renormalize', action='store_true',
		help='Renormalize abundances to 100 percent at each rank,\
				for instance if an organism has a species but not genus label.')
	parser.add_argument('--output', default='abundances.tsv',
		help='Output abundances file. Default: abundances.tsv')
	parser.add_argument('--pct_id', type=float, default=0.5,
		help='Minimum percent identity from reference to count a hit.')
	parser.add_argument('--no_quantify_unmapped', action='store_true',
		help='Do not factor in unmapped reads in abundance estimation.')
	parser.add_argument('--read_cutoff', type=int, default=1,
		help='Number of reads to count an organism as present.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to input file name(s).')
	parser.add_argument('--strain_level', action='store_true',
		help='Use this flag to profile strains (off by default).')
	parser.add_argument('--temp_dir', default = 'TEMP_metalign/',
		help='Directory to write temporary files to.')
	parser.add_argument('--verbose', action='store_true',
		help='Print verbose output.')
	args = parser.parse_args()
	return args


def main():
	args = metalign_parseargs()
	if not args.temp_dir.endswith('/'):
		args.temp_dir += '/'
	if not args.data.endswith('/'):
		args.data += '/'
	# Set arguments that default to AUTO
	if args.dbinfo_in == 'AUTO':
		args.data + 'db_info.txt'
	if args.db_dir == 'AUTO':
		args.data + 'organism_files/'

	# Ensure that arguments agree between scripts
	args.db = args.temp_dir + 'cmashed_db.fna'
	args.dbinfo = args.temp_dir + 'subset_db_info.txt'
	args.dbinfo_out = args.dbinfo
	args.infiles = [args.reads]  # map_and_profile expects a list
	args.cmash_results = 'NONE'

	# Run the database selection and map/profile routines
	select.select_main(args)  # runs select_db routine
	mapper.map_main(args)  # runs map_and_profile routine


if __name__ == '__main__':
	main()
#
