import argparse, os, shlex, subprocess, sys, time


# Globals: Program timer, download links, and CAMI-relevant taxonomic ranks
start = time.time()
ASM_ORDER = ['archaea', 'bacteria', 'fungi', 'protozoa', 'viral']
DB_ORDER = ['refseq', 'genbank']
TAXDUMP_LOC = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip'
# Links relevant CAMI ranks to an index in a lineage list (see trace_lineages)
CAMI_RANKS = {'superkingdom':0, 'phylum':1, 'class':2, 'order':3,
				'family':4, 'genus':5, 'species':6, 'strain':7}


def echo(msg, prepend='', postpend=''):
        global start
        seconds = time.time() - start
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        hms = "%02d:%02d:%02d" % (h, m, s)
        print(prepend + '['+hms+'] ' + msg + postpend)


# Credit (Original code under MIT License):
# https://stackoverflow.com/questions/6169217/replace-console-output-in-python
def progressBar(value, endvalue, bar_length=50):
	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rProgress: [{0}] {1}%".format(arrow + spaces,
		int(round(percent * 100))))
	sys.stdout.flush()


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(
		description='Download all RefSeq genomes and extract needed info.')
	parser.add_argument('--dir', default='micopdb/',
		help = 'Name of directory to store database and info.')
	parser.add_argument('--mode', default='all',
		choices=['all', 'setup', 'dl-genomes', 'build_info'],
		help = 'What steps of database building to run. Default: all')
	parser.add_argument('--no_clean', action='store_true',
		help = 'Do not clean up intermediate files after running.')
	parser.add_argument('--no_combine', action='store_true',
		help = 'Do not combine individual organism files into one database.')
	parser.add_argument('--test', action='store_true',
		help = 'Do a small test run (should take about 1-3 minutes).')
	parser.add_argument('--verbose', action='store_true',
		help = 'Print progress updates (otherwise, prints nothing).')
	args = parser.parse_args()
	return args


def download_summaries(args):
	spec_fnames, url_prefix = [], 'ftp://ftp.ncbi.nlm.nih.gov/genomes/'
	for i in range(10):
		# Iterate through clades and refseq/genbank, create specific filename
		dbtype, clade = DB_ORDER[int(i/5)], ASM_ORDER[i % 5]
		link = url_prefix + dbtype + '/' + clade + '/assembly_summary.txt'
		spec_fnames.append(dbtype + '_' + clade + '_' + 'assembly_summary.txt')
		with(open('logfile', 'a')) as log:
			subprocess.Popen(['wget', link], stdout=log, stderr=log).wait()

		# If this is a test run, retain only 10 lines per clade
		if args.test:
			with(open('TEMP', 'w')) as tmp:
				subprocess.Popen(['head', 'assembly_summary.txt'],
					stdout=tmp).wait()
				subprocess.Popen(['rm', 'assembly_summary.txt']).wait()
				subprocess.Popen(['mv', 'TEMP', 'assembly_summary.txt']).wait()

		# Move assembly_summary.txt to a unique filename to avoid overwriting
		subprocess.Popen(['mv', 'assembly_summary.txt', spec_fnames[-1]]).wait()

	# Combine all refseq entries and unique genbank entries into one single file
	with(open('assembly_summary.txt', 'w')) as outfile:
		for specfile in spec_fnames:
			rfsq = 'refseq' in specfile
			acc, prev_acc = '', ''  # used to check duplicate entries
			with(open(specfile, 'r')) as infile:
				infile.readline(); infile.readline()  # skip header lines
				for line in infile:  #write all refseq, unpaired genbank entries
					if rfsq or not ('identical' in line or 'different' in line):
						acc = line.split('\t')[0]
						if acc != prev_acc:
							outfile.write(line)
						prev_acc = acc


def setup(args):
	download_summaries(args)  # download and setup assembly_summary.txt file
	# Now extract relevant file names to download from the assembly summaries
	cmd = 'awk -F "\t" \'$11=="latest"{print $20}\' assembly_summary.txt'
	with(open('ftpdirpaths', 'w')) as cmdfile:
		subprocess.Popen(shlex.split(cmd), stdout = cmdfile).wait()
	cmd = 'awk \'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}\
		{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}\' \
		ftpdirpaths'
	with(open('ftpfilepaths', 'w')) as cmdfile:
		subprocess.Popen(shlex.split(cmd), stdout = cmdfile).wait()
	cmd = 'awk \'BEGIN{FS=OFS="/";filesuffix="assembly_report.txt"}\
		{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}\' \
		ftpdirpaths'
	with(open('ftpreportpaths', 'w')) as cmdfile:
		subprocess.Popen(shlex.split(cmd), stdout = cmdfile).wait()

	# Download and extract taxonomy files
	with(open('logfile', 'a')) as logfile:
		subprocess.Popen(['wget', TAXDUMP_LOC],
			stdout=logfile, stderr=logfile).wait()
		if not os.path.isdir('taxonomy/'):
			os.mkdir('taxonomy/')
		subprocess.Popen(['mv', 'new_taxdump.zip', 'taxonomy/']).wait()
		os.chdir('taxonomy/')
		subprocess.Popen(['unzip', 'new_taxdump.zip'],
			stdout=logfile, stderr=logfile).wait()
		os.chdir('..')
		subprocess.Popen(['cp', 'taxonomy/names.dmp', '.']).wait()
		subprocess.Popen(['cp', 'taxonomy/nodes.dmp', '.']).wait()

	# Finally, clean some old files out unless --no_clean is used
	if not args.no_clean:
		for fn in spec_fnames:
			subprocess.Popen(['rm', fn]).wait()
		subprocess.Popen(['rm', 'assembly_summary.txt', 'ftpdirpaths']).wait()
		subprocess.Popen(['rm', '-r', 'taxonomy/']).wait()


def build_taxtree(args):
	# taxtree is stored as a dict with TaxID keys linking to
	# 		taxonomic name, taxonomic rank, and parent TaxID
	# Thus, it will be easy to traverse through lineage given starting TaxID
	taxtree = {}

	with(open('names.dmp', 'r')) as names:
		for line in names:
			if 'scientific name' not in line:
				continue  # skip alternate names and other unnecessary info
			#splits = line.split('|')[1].strip()
			taxid = line.split()[0]
			name = line.split('|')[1].strip()
			taxtree[taxid] = [name]  # Key is TaxID, value is tax. name

	with(open('nodes.dmp', 'r')) as nodes:
		for line in nodes:
			splits = line.split()
			# add taxonomic rank and parent TaxID to list of info for this taxID
			taxtree[splits[0]].extend([splits[4], splits[2]])
	return taxtree


def download_genomes(args):
	# clear files because they will be appended to
	open('ref_db.fasta', 'w').close()
	open('reports.txt', 'w').close()

	# number of downloads needed and how many are done, for progress tracking
	done, num_dl = 0, int(subprocess.check_output(['wc', '-l',
		'ftpfilepaths']).split()[0])

	# We now open the file handlers for appending, download the files,
	# 		and append them to the growing reference database and info file
	if not args.no_combine:
		refseq = open('ref_db.fasta', 'a')
	nc = '../' if args.no_combine else ''  # add to filenames if args.no_combine
	reports = open(nc+'reports.txt', 'a')  # temp file, db_info constructed later
	if args.no_combine:  # if keeping individual org. files, mkdir for them
		os.makedirs('organism_files')
		os.chdir('organism_files')

	with(open(nc+'ftpfilepaths', 'r')) as ftpfilepaths:
		with(open(nc+'ftpreportpaths', 'r')) as ftpreportpaths:
			for fileline in ftpfilepaths:
				reportline = ftpreportpaths.readline()
				# download the files
				with(open(nc+'logfile', 'a')) as logfile:
					subprocess.Popen(['wget', fileline.strip()],
						stdout=logfile, stderr=logfile).wait()
					subprocess.Popen(['wget', reportline.strip()],
						stdout=logfile, stderr=logfile).wait()

				if not args.no_combine:
					# grab the file names, decompress the genome, and append to
					# 		ref_db.fasta and dbinfo.txt
					fname = fileline.strip().split('/')[-1]
					p = subprocess.Popen(['cat', fname], stdout=subprocess.PIPE)
					subprocess.Popen(['gunzip'],
						stdin=p.stdout, stdout=refseq).wait()
					subprocess.Popen(['rm', fname]).wait()
				# always combine report files
				reportname = reportline.strip().split('/')[-1]
				subprocess.Popen(['cat', reportname], stdout=reports).wait()
				subprocess.Popen(['rm', reportname]).wait()
				if args.verbose:
					done += 1
					progressBar(done, num_dl)
	reports.close()
	if args.no_combine:
		os.chdir('..')
	else:
		refseq.close()


def get_accession_lengths(args):
	# reads through ref_db.fasta, matching accession to entry length in bases
	acc2len = {}
	accession, curlen = 'IGNORE', 0
	with(open('ref_db.fasta', 'r')) as refseq:
		for line in refseq:
			if line.startswith('>'):
				acc2len[accession] = str(curlen)
				accession = line.split()[0][1:]
				curlen = 0
			else:
				curlen += len(line.strip())
	del acc2len['IGNORE']  # remove initial entry created by above process
	return acc2len


def trace_lineages(taxid, taxtree):
	# Uses taxtree to find lineage of taxid in both taxonomic names and IDs
	name_lineage, taxid_lineage = ['' for i in range(8)], ['' for i in range(8)]
	cur_taxid = taxid

	# Record lowest level taxonomic info if it's below species level
	if cur_taxid in taxtree:
		name, rank, parent = taxtree[cur_taxid]
	else:
		return 'NONE', 'NONE'
	if rank not in CAMI_RANKS:  # strains are recorded as 'no rank' in nodes.dmp
		name_lineage[-1] = name
		taxid_lineage[-1] = cur_taxid
		cur_taxid = parent  # traverse up to parent taxid

	# Now traverse up the tree
	while cur_taxid != '1':  # while we're not at the root
		if cur_taxid in taxtree:
			name, rank, parent = taxtree[cur_taxid]
		else:
			return 'NONE', 'NONE'
		if rank in CAMI_RANKS:  # if this is a relevant CAMI taxonomic rank
			index = CAMI_RANKS[rank]  # then record this node's info at the rank
			name_lineage[index] = name
			taxid_lineage[index] = cur_taxid
		cur_taxid = parent  # traverse up to parent taxid and repeat
	return '|'.join(name_lineage), '|'.join(taxid_lineage)


def construct_dbinfo(args, taxtree):
	# Here we parse the raw reports in reports.txt, extract the
	# 		TaxID and Accessions, and reconstruct the lineage for the TaxID
	# 		using the taxtree. Then we store all of this info in db_info.txt
	acc2len = get_accession_lengths(args)  # dict matching accession to length
	taxid, accessions, name_lin, taxid_lin = '', [], '', ''

	# number of downloads needed and how many are done, for progress tracking
	done, num_dl = 0, int(subprocess.check_output(['grep', '-c',
		'^# Taxid', 'reports.txt']).strip())

	with(open('reports.txt', 'r')) as reports:
		with(open('db_info.txt', 'w')) as dbinfo:
			dbinfo.write('Accesion\tLength   \tTaxID   \t')
			dbinfo.write('Lineage \tTaxID_Lineage\n')
			for line in reports:
				if line.startswith('# Assembly name:'):  # new report
					for acc in accessions:  # write info to db_info for each acc
						if name_lin != 'NONE' and acc in acc2len:
							info=[acc, acc2len[acc], taxid, name_lin, taxid_lin]
							dbinfo.write('\t'.join(info) + '\n')
					taxid, accessions = '', []
					if args.verbose:
						done += 1
						progressBar(done, num_dl)

				elif line.startswith('# Taxid:'):  # grab taxid, trace lineage
					taxid = line.strip().split()[-1]
					name_lin, taxid_lin = trace_lineages(taxid, taxtree)

				elif not line.startswith('#'):  # line with an accession number
					accessions.append(line.split()[6])
				# All other lines are unnecessary


def remove_dupl_dbinfo(args):  # removes duplicate entries in db_info
	prev = ''  # keeps track of previous accession; if same as current, remove
	with(open('db_info.txt', 'r')) as infile:
		with(open('db_info-clean.txt', 'w')) as outfile:
			for line in infile:
				cur = line.split()[0]
				if cur != prev:
					outfile.write(line)
				prev = cur
	subprocess.Popen(['rm', 'db_info.txt']).wait()
	subprocess.Popen(['mv', 'db_info-clean.txt', 'db_info.txt']).wait()


def clean_refseq(args):
	info_accessions = {}  # stores accessions from db_info
	write = True  # whether to write out current accession info

	with(open('db_info.txt', 'r')) as dbinfo:
		dbinfo.next()
		for line in dbinfo:
			acc = line.split()[0]
			info_accessions[acc] = 0

	with(open('ref_db.fasta', 'r')) as infile:
		with(open('refseq-clean.fasta', 'w')) as outfile:
			for line in infile:
				if line.startswith('>'):
					acc = line[1:].split()[0]
					write = (acc in info_accessions)
				if write:
					outfile.write(line)
	subprocess.Popen(['rm', 'ref_db.fasta']).wait()
	subprocess.Popen(['mv', 'refseq-clean.fasta', 'ref_db.fasta']).wait()


def main():
	# Parse user arguments, set up working directory and logfile
	args = parseargs()
	if not args.dir.endswith('/'):
		args.dir += '/'
	if not os.path.isdir(args.dir):
		os.makedirs(args.dir)
	os.chdir(args.dir)
	open('logfile', 'w').close()  # clear the logfile

	# Call the methods to construct ref_db.fasta and db_info.txt
	if args.mode == 'all' or args.mode == 'setup':
		if args.verbose:
			echo('Downloading filepaths and taxonomy dump (<5 mins)...')
		setup(args)
	if args.mode != 'setup':
		if args.verbose:
			echo('Building taxonomic tree (<1 min)...')
		taxtree = build_taxtree(args)
		if args.mode == 'all' or args.mode == 'dl-genomes':
			if args.verbose:
				echo('Downloading genomes and building database (very long)...')
			download_genomes(args)
		if (args.mode == 'all' or args.mode == 'build_info'):
			if not args.no_combine:
				if args.verbose:
					echo('Building database info file (~1 hour)...',
						prepend='\n')
				construct_dbinfo(args, taxtree)
				if args.verbose:
					echo('Adjusting databases (<1 hour)...', prepend='\n')
				remove_dupl_dbinfo(args)
				clean_refseq(args)

	if not args.no_clean:
		# Clear temporary files
		subprocess.Popen(['rm', 'ftpfilepaths', 'ftpreportpaths',
		 	'reports.txt', 'logfile']).wait()
		subprocess.Popen(['rm', 'names.dmp', 'nodes.dmp']).wait()


if __name__ == '__main__':
	main()
#
