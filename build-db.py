import argparse, os, shlex, subprocess, sys, time


# Globals: Program timer, download links, and CAMI-relevant taxonomic ranks
start = time.time()
ASM = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt']
ASM_ORDER = ['archaea', 'bacteria', 'fungi', 'protozoa', 'viral']
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
	parser.add_argument('--test', action='store_true',
		help = 'Do a small test run (should take about 1-3 minutes).')
	parser.add_argument('--verbose', action='store_true',
		help = 'Print progress updates (otherwise, prints nothing).')
	args = parser.parse_args()
	return args


def download_info(args):
	clade_fnames = [o + '_' + 'assembly_summary.txt' for o in ASM_ORDER]
	for i in range(len(ASM)):
		# Download assembly summary file with genome ftp addresses
		with(open('logfile', 'a')) as logfile:
			subprocess.Popen(['wget', ASM[i]],
				stdout=logfile, stderr=logfile).wait()

		# If this is a test run, retain only 10 lines per clade
		if args.test:
			with(open('tmpfile', 'w')) as tmpfile:
				subprocess.Popen(['head', 'assembly_summary.txt'],
					stdout=tmpfile).wait()
				subprocess.Popen(['rm', 'assembly_summary.txt']).wait()
				subprocess.Popen(['mv', 'tmpfile',
					'assembly_summary.txt']).wait()

		# Move assembly_summary.txt to a clade-specific file name so it's not
		# 		overwritten when we download the next one
		subprocess.Popen(['mv', 'assembly_summary.txt', clade_fnames[i]]).wait()

	# Now we cat these clade-specific summaries into one single file
	cmd = ['cat'] + clade_fnames
	with(open('assembly_summary.txt', 'w')) as cmdfile:
		subprocess.Popen(cmd, stdout = cmdfile).wait()

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

	# Finally, clean some old files out
	for fn in clade_fnames:
		subprocess.Popen(['rm', fn]).wait()
	subprocess.Popen(['rm', 'assembly_summary.txt', 'ftpdirpaths']).wait()
	subprocess.Popen(['cp', 'taxonomy/names.dmp', '.']).wait()
	subprocess.Popen(['cp', 'taxonomy/nodes.dmp', '.']).wait()
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
	open('refseq.fasta', 'w').close()
	open('reports.txt', 'w').close()

	# number of downloads needed and how many are done, for progress tracking
	done, num_dl = 0, int(subprocess.check_output(['wc', '-l',
		'ftpfilepaths']).split()[0])

	# We now open the file handlers for appending, download the files,
	# 		and append them to the growing reference database and info file
	refseq = open('refseq.fasta', 'a')
	reports = open('reports.txt', 'a')  # temp file, db_info constructed later
	with(open('ftpfilepaths', 'r')) as ftpfilepaths:
		with(open('ftpreportpaths', 'r')) as ftpreportpaths:
			for fileline in ftpfilepaths:
				reportline = ftpreportpaths.readline()
				# download the files
				with(open('logfile', 'a')) as logfile:
					subprocess.Popen(['wget', fileline.strip()],
						stdout=logfile, stderr=logfile).wait()
					subprocess.Popen(['wget', reportline.strip()],
						stdout=logfile, stderr=logfile).wait()

				# grab the file names, decompress the genome, and append to
				# 		refseq.fasta and dbinfo.txt
				filename = fileline.strip().split('/')[-1]
				reportname = reportline.strip().split('/')[-1]
				p = subprocess.Popen(['cat', filename], stdout=subprocess.PIPE)
				subprocess.Popen(['gunzip'],
					stdin=p.stdout, stdout=refseq).wait()
				subprocess.Popen(['cat', reportname], stdout=reports).wait()

				# Finally, remove the individual files and update progress
				subprocess.Popen(['rm', filename, reportname]).wait()
				if args.verbose:
					done += 1
					progressBar(done, num_dl)
	refseq.close()
	reports.close()


def get_accession_lengths(args):
	# reads through refseq.fasta, matching accession to entry length in bases
	acc2len = {}
	accession, curlen = 'IGNORE', 0
	with(open('refseq.fasta', 'r')) as refseq:
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

	with(open('refseq.fasta', 'r')) as infile:
		with(open('refseq-clean.fasta', 'w')) as outfile:
			for line in infile:
				if line.startswith('>'):
					acc = line[1:].split()[0]
					write = (acc in info_accessions)
				if write:
					outfile.write(line)
	subprocess.Popen(['rm', 'refseq.fasta']).wait()
	subprocess.Popen(['mv', 'refseq-clean.fasta', 'refseq.fasta']).wait()


def main():
	# Parse user arguments, set up working directory and logfile
	args = parseargs()
	if not args.dir.endswith('/'):
		args.dir += '/'
	if not os.path.isdir(args.dir):
		os.makedirs(args.dir)
	os.chdir(args.dir)
	open('logfile', 'w').close()  # clear the logfile

	# Call the methods to construct refseq.fasta and db_info.txt
	if args.verbose:
		echo('Downloading filepaths and taxonomy dump (<5 mins)...')
	download_info(args)
	if args.verbose:
		echo('Building taxonomic tree (<1 min)...')
	taxtree = build_taxtree(args)
	if args.verbose:
		echo('Downloading genomes and building database (very long)...')
	download_genomes(args)
	if args.verbose:
		echo('Building database info file (~1 hour)...', prepend='\n')
	construct_dbinfo(args, taxtree)
	if args.verbose:
		echo('Cleaning files (<1 hour)...', prepend='\n')
	remove_dupl_dbinfo(args)
	clean_refseq(args)

	# Clear temporary files
	subprocess.Popen(['rm', 'ftpfilepaths', 'ftpreportpaths',
	 	'reports.txt', 'logfile']).wait()
	subprocess.Popen(['rm', 'names.dmp', 'nodes.dmp']).wait()


if __name__ == '__main__':
	main()
#

