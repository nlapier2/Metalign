import argparse, math, os, subprocess, sys


__location__ = os.path.realpath(os.path.join(os.getcwd(),
								os.path.dirname(__file__))) + '/'


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Run CMash and" +
				" select a subset of the whole database to align to.")
	parser.add_argument('reads', help='Path to reads file.')
	parser.add_argument('db_dir',
		help='Directory with all organism files in the full database.')
	parser.add_argument('kmc_dir',
		help='Directory to write temporary KMC files to.')
	parser.add_argument('--cmash_results', default='NONE',
		help='Can specfily location of CMash query results if already done.')
	parser.add_argument('--cutoff', type=float, default=-1.0,
		help='CMash cutoff value. Default is 1/(log10(reads file bytes)**2).')
	parser.add_argument('--dbinfo_in', default='AUTO',
		help='Specify location of db_info file. Default is data/db_info.txt')
	parser.add_argument('--dbinfo_out', default='AUTO',
		help='Where to write subset db_info. Default: data/subset_db_info.txt')
	parser.add_argument('--keep_kmc_files', action = 'store_true',
		help='Keep KMC files instead of deleting after this script finisishes.')
	parser.add_argument('--output', default='cmashed_db.fna',
		help='Path to where to write the output database.')
	parser.add_argument('--strain_level', action='store_true',
		help='Include all strains above cutoff. Default: 1 strain per species.')
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
	kmc_loc = __location__ + 'kmc/kmc'
	db_60mers_loc = __location__ + 'data/cmash_db_n1000_k60_dump'

	subprocess.Popen([kmc_loc, '-v', '-k60', '-m200', '-sm', '-fq', '-ci2',
		'-cs3', '-t48', '-jlog_sample', args.reads,
		args.kmc_dir + 'reads_60mers', args.kmc_dir]).wait()

	subprocess.Popen([kmc_loc+'_tools', 'simple', db_60mers_loc,
		args.kmc_dir + 'reads_60mers', 'intersect',
		args.kmc_dir + '60mers_intersection']).wait()

	subprocess.Popen([kmc_loc+'_dump', args.kmc_dir + '60mers_intersection',
		args.kmc_dir + '60mers_intersection_dump']).wait()

	with(open(args.kmc_dir + '60mers_intersection_dump', 'r')) as infile:
		with(open(args.kmc_dir + '60mers_intersection_dump.fa', 'w')) as fasta:
			for line in infile:
				seq = line.split()[0]
				fasta.write('>seq' + '\n' + seq + '\n')


	#cat_proc = subprocess.Popen(['cat', args.kmc_dir+'60mers_intersection_dump'],
	#	stdout=subprocess.PIPE)
	#cut_proc = subprocess.Popen(['cut', '-f1'], stdin=cat_proc.stdout,
	#	stdout=subprocess.PIPE)
	#with(open(args.kmc_dir + '60mers_intersection_dump.fa', 'w')) as kmc_res:
	#	sed_proc = subprocess.Popen(['sed', "'s/^/>seq\n/g'"],
	#		stdin = cut_proc.stdout, stdout = kmc_res).wait()


def run_cmash_and_cutoff(args, taxid2info):
	cmash_db_loc = __location__ + 'data/cmash_db_n1000_k60.h5'
	cmash_filter_loc = __location__ + 'data/cmash_filter_n1000_k60_30-60-10.bf'
	if args.cmash_results == 'NONE':
		# temporarily write cmash results to --output file; should be writeable
		cmash_out = args.output
		script_loc = __location__ + 'CMash/scripts/StreamingQueryDNADatabase.py'
		cmash_proc = subprocess.Popen(['python', script_loc,
			args.kmc_dir + '60mers_intersection_dump.fa', cmash_db_loc,
			args.output, '30-60-10', '-c', '0', '-r', '1000000', '-v',
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
	open(args.output, 'w').close()  # clear cmash results; no longer needed
	with(open(args.output, 'a')) as outfile:
		for organism in organisms_to_include:
			organism_fname = args.db_dir + organism
			# write organisms to full db via cat to append-mode file handler
			subprocess.call(['cat', organism_fname], stdout=outfile)

	with(open(args.dbinfo_out, 'w')) as outfile:
		# write header lines
		outfile.write('Accesion\tLength\tTaxID\tLineage\tTaxID_Lineage\n')
		outfile.write('Unmapped\t0\tUnmapped\t|||||||Unmapped\t|||||||Unmapped\n')
		for organism in organisms_to_include:
			taxid = organism.split('taxid_')[1].split(
				'_genomic.fna')[0].replace('_', '.')
			len = taxid2info[taxid][1]
			namelin, taxlin = taxid2info[taxid][2], taxid2info[taxid][3]
			for acc in taxid2info[taxid][0]:
				outfile.write('\t'.join([acc,len,taxid,namelin,taxlin]) + '\n')


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
	if not args.kmc_dir.endswith('/'):
		args.kmc_dir += '/'
	if args.dbinfo_in == 'AUTO':
		args.dbinfo_in = __location__ + 'data/db_info.txt'
	if args.dbinfo_out == 'AUTO':
		args.dbinfo_out = __location__ + 'data/subset_db_info.txt'

	taxid2info = read_dbinfo(args)
	if args.cmash_results == 'NONE':
		run_kmc_steps(args)
	organisms_to_include = run_cmash_and_cutoff(args, taxid2info)
	make_db_and_dbinfo(args, organisms_to_include, taxid2info)

	if not args.keep_kmc_files:
		subprocess.Popen(['rm', args.kmc_dir + 'reads_60mers.kmc_pre',
		args.kmc_dir + 'reads_60mers.kmc_suf',
		args.kmc_dir + '60mers_intersection.kmc_pre',
		args.kmc_dir + '60mers_intersection.kmc_suf',
		args.kmc_dir + '60mers_intersection_dump',
		args.kmc_dir + '60mers_intersection_dump.fa']).wait()


if __name__ == '__main__':
	main()
#
