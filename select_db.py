import argparse, math, os, subprocess, sys, tempfile


__location__ = os.path.realpath(os.path.join(os.getcwd(),
								os.path.dirname(__file__))) + '/'


def select_parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Run CMash and" +
				" select a subset of the whole database to align to.")
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('data',
		help='Path to data/ directory with the files from setup_data.sh')
	parser.add_argument('--cmash_results', default='NONE',
		help='Can specify location of CMash query results if already done.')
	parser.add_argument('--cutoff', type=float, default=-1.0,
		help='CMash cutoff value. Default is 1/(log10(reads file bytes)**2).')
	parser.add_argument('--db', default='cmashed_db.fna',
		help='Path to where to write the output database.')
	parser.add_argument('--db_dir', default='AUTO',
		help='Directory with all organism files in the full database.')
	parser.add_argument('--dbinfo_in', default='AUTO',
		help='Specify location of db_info file. Default is data/db_info.txt')
	parser.add_argument('--dbinfo_out', default='AUTO',
		help='Where to write subset db_info. Default: data/subset_db_info.txt')
	parser.add_argument('--input_type', default='AUTO',
		choices=['fastq', 'fasta', 'AUTO'],
		help='Type of input file (fastq/fasta). Default: try to auto-determine')
	parser.add_argument('--keep_temp_files', action = 'store_true',
		help='Keep KMC files instead of deleting after this script finishes.')
	parser.add_argument('--strain_level', action='store_true',
		help='Include all strains above cutoff. Default: 1 strain per species.')
	parser.add_argument('--temp_dir', default = 'AUTO/',
		help='Directory to write temporary files to.')
	parser.add_argument('--threads', type=int, default=4,
		help='How many compute threads for KMC to use. Default: 4')
	args = parser.parse_args()
	return args


def read_dbinfo(args):
	taxid2info = {}
	with(open(args.dbinfo_in, 'r')) as infile:
		infile.readline()  # skip header line
		for line in infile:
			splits = line.strip().split('\t')
			acc, taxid = splits[0], splits[2]
			if taxid not in taxid2info:
				# first element stores all accessions for this taxid
				taxid2info[taxid] = [[splits[0]], splits[1]]
				taxid2info[taxid].extend(splits[3:])
			else:
				taxid2info[taxid][0].append(acc)
	return taxid2info


def run_kmc_steps(args):
	kmc_loc = __location__ + 'KMC/bin/kmc'
	db_60mers_loc = args.data + 'cmash_db_n1000_k60_dump'
	if args.input_type == 'fastq':
		type_arg = '-fq'
	else:
		type_arg = '-fa'

	subprocess.Popen([kmc_loc, '-v', '-k60', type_arg, '-ci2', '-cs3',
		'-t' + str(args.threads), '-jlog_sample', args.reads,
		args.temp_dir + 'reads_60mers', args.temp_dir]).wait()

	subprocess.Popen([kmc_loc+'_tools', 'simple', db_60mers_loc,
		args.temp_dir + 'reads_60mers', 'intersect',
		args.temp_dir + '60mers_intersection']).wait()

	subprocess.Popen([kmc_loc+'_dump', args.temp_dir + '60mers_intersection',
		args.temp_dir + '60mers_intersection_dump']).wait()

	with(open(args.temp_dir + '60mers_intersection_dump', 'r')) as infile:
		with(open(args.temp_dir + '60mers_intersection_dump.fa', 'w')) as fasta:
			for line in infile:
				seq = line.split()[0]
				fasta.write('>seq' + '\n' + seq + '\n')


def run_cmash_and_cutoff(args, taxid2info):
	cmash_db_loc = args.data + 'cmash_db_n1000_k60.h5'
	cmash_filter_loc = args.data + 'cmash_filter_n1000_k60_30-60-10.bf'
	if args.cmash_results == 'NONE':
		cmash_out = args.temp_dir + 'cmash_query_results.csv'
		script_loc = __location__ + 'CMash/scripts/StreamingQueryDNADatabase.py'
		cmash_proc = subprocess.Popen(['python', script_loc,
			args.temp_dir + '60mers_intersection_dump.fa', cmash_db_loc,
			cmash_out, '30-60-10', '-c', '0', '-r', '1000000', '-v',
			'-f', cmash_filter_loc, '--sensitive']).wait()
	else:
		cmash_out = args.cmash_results

	organisms_to_include, species_included = [], {}
	with(open(cmash_out, 'r')) as cmash_results:
		cmash_results.readline()  # skip header line
		for line in cmash_results:
			splits = line.strip().split(',')
			organism, containment_index = splits[0], float(splits[-1])
			if containment_index >= args.cutoff:
				if not args.strain_level:
					taxid = organism.split('taxid_')[1].split(
						'_genomic.fna')[0].replace('_', '.')
					species = taxid2info[taxid][3].split('|')[-2]
					if species not in species_included or species == '':
						species_included[species] = 1
					else:
						continue
				organisms_to_include.append(organism)
	return organisms_to_include


def make_db_and_dbinfo(args, organisms_to_include, taxid2info):
	open(args.db, 'w').close()  # clear cmash results; no longer needed
	with(open(args.db, 'a')) as outfile:
		for organism in organisms_to_include:
			organism_fname = args.db_dir + organism
			# write organisms to full db via cat to append-mode file handler
			subprocess.Popen(['zcat', organism_fname], stdout=outfile).wait()

	with(open(args.dbinfo_out, 'w')) as outfile:
		# write header lines
		outfile.write('Accesion\tLength\tTaxID\tLineage\tTaxID_Lineage\n')
		outfile.write('Unmapped\t0\tUnmapped\t|||||||Unmapped\t|||||||Unmapped\n')
		for organism in organisms_to_include:
			taxid = organism.split('taxid_')[1].split(
				'_genomic.fna')[0].replace('_', '.')
			length = taxid2info[taxid][1]
			namelin, taxlin = taxid2info[taxid][2], taxid2info[taxid][3]
			for acc in taxid2info[taxid][0]:
				outfile.write('\t'.join([acc,length,taxid,namelin,taxlin]) + '\n')


def select_main(args = None):
	if args == None:
		args = select_parseargs()
	if args.cutoff == -1.0:  # not set by user
		fsize_bytes = os.path.getsize(args.reads)
		args.cutoff = 1.0 / (math.log10(fsize_bytes) ** 2)
	elif args.cutoff < 0.0 or args.cutoff > 1.0:
		print('Error: args.cutoff must be between 0 and 1, inclusive.')
		sys.exit()
	if not args.data.endswith('/'):
		args.data += '/'
	if args.db_dir == 'AUTO':
		args.db_dir = args.data + 'organism_files/'
	if not args.db_dir.endswith('/'):
		args.db_dir += '/'
	if args.temp_dir == 'AUTO/':
		args.temp_dir = tempfile.mkdtemp(prefix=__location__+'data/')
	if not args.temp_dir.endswith('/'):
		args.temp_dir += '/'
	if not os.path.exists(args.temp_dir):
		os.makedirs(args.temp_dir)
	if args.dbinfo_in == 'AUTO':
		args.dbinfo_in = args.data + 'db_info.txt'
	if args.dbinfo_out == 'AUTO':
		args.dbinfo_out = args.data + 'subset_db_info.txt'
	if args.input_type == 'AUTO':
		splits = args.reads.split('.')
		if splits[-1] == 'gz':  # gz doesn't help determine file type
			splits = splits[:-1]
		if splits[-1] in ['fq', 'fastq']:
			args.input_type = 'fastq'
		elif splits[-1] in ['fa', 'fna', 'fasta']:
			args.input_type = 'fasta'
		else:
			sys.exit('Could not auto-determine file type. Use --input_type.')

	taxid2info = read_dbinfo(args)
	if args.cmash_results == 'NONE':
		run_kmc_steps(args)
	organisms_to_include = run_cmash_and_cutoff(args, taxid2info)
	make_db_and_dbinfo(args, organisms_to_include, taxid2info)

	if not args.keep_temp_files and args.cmash_results == 'NONE':
		subprocess.Popen(['rm', args.temp_dir + 'reads_60mers.kmc_pre',
		args.temp_dir + 'reads_60mers.kmc_suf',
		args.temp_dir + '60mers_intersection.kmc_pre',
		args.temp_dir + '60mers_intersection.kmc_suf',
		args.temp_dir + '60mers_intersection_dump',
		args.temp_dir + '60mers_intersection_dump.fa']).wait()


if __name__ == '__main__':
	args = select_parseargs()
	select_main(args)
#
