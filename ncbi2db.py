import argparse, glob, gzip, os, subprocess


# define globals
TAXDUMP_LOC = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip'
CAMI_RANKS = {'superkingdom':0, 'phylum':1, 'class':2, 'order':3,
				'family':4, 'genus':5, 'species':6, 'strain':7}


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given genome &" +
				" assembly report files from NCBI, build MiCoP2 database.")
	parser.add_argument('--input_dir',default='ncbi_rsync_all_genomes_reports/',
		help = 'Name of directory with genome/report files from NCBI.')
	parser.add_argument('--output_dir', default='organism_files/',
		help = 'Name of directory to write selected genome files to.')
	args = parser.parse_args()
	return args


def build_taxtree():
	# taxtree is stored as a dict with TaxID keys linking to
	# 		taxonomic name, taxonomic rank, and parent TaxID
	# Thus, it will be easy to traverse through lineage given starting TaxID
	taxtree = {}

	with(open('taxonomy/names.dmp', 'r')) as names:
		for line in names:
			if 'scientific name' not in line:
				continue  # skip alternate names and other unnecessary info
			#splits = line.split('|')[1].strip()
			taxid = line.split()[0]
			name = line.split('|')[1].strip()
			taxtree[taxid] = [name]  # Key is TaxID, value is tax. name

	with(open('taxonomy/nodes.dmp', 'r')) as nodes:
		for line in nodes:
			splits = line.split()
			# add taxonomic rank and parent TaxID to list of info for this taxID
			taxtree[splits[0]].extend([splits[4], splits[2]])
	return taxtree


# Uses taxtree to find lineage of taxid in both taxonomic names and IDs
def trace_lineages(taxid, taxtree):
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


def get_taxonomy_info():
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

		# build taxid/name lineage tree, and taxid to full name lineage
		taxtree = build_taxtree()
		taxid2namelin = {}
		with(open('taxonomy/fullnamelineage.dmp', 'r')) as infile:
			for line in infile:
				splits = line.strip().split('|')
				#taxid, lineage = splits[0].strip(), splits[-2].split('; ')
				#taxid2namelin[taxid] = lineage
				taxid2namelin[splits[0].strip()] = splits[-2]
		subprocess.Popen(['rm', '-r', 'taxonomy/']).wait()  # clean tax. files
		return taxtree, taxid2namelin


# Searches assembly reports in input_dir, checks to see if microbial,
#  	and if they are, puts in mapping of taxid to all matching organisms.
#  	Also keeps track of all versions of an organism and final one to keep.
def get_taxids2asmnames(args, taxid2namelin):
	ignored_taxa = ['Metazoa', 'Embryophyta', 'unclassified sequences']
	taxid2asmnames, name2final_ver = {}, {}
	# only check reports that have a matching genomic file
	for genomic_name in glob.glob(args.input_dir+'*_genomic.fna.gz'):
		asm_acc = genomic_name.split('/')[-1].split('_genomic.fna.gz')[0]
		org_name, taxid = asm_acc[4:13], ''
		with(open(args.input_dir + asm_acc +'_assembly_report.txt',
					'r')) as report_file:
			for line in report_file:
				if 'Taxid' in line:
					taxid = line.strip().split()[-1]
					break  # we have the info we need

		# filter out taxid if non-microbial (animal/plant/unassigned)
		if taxid not in taxid2namelin:
			continue
		lin = taxid2namelin[taxid]
		if any([i in lin for i in ignored_taxa]):
			continue
		if taxid in taxid2asmnames:
			taxid2asmnames[taxid].append([asm_acc, org_name])
		else:
			taxid2asmnames[taxid] = [[asm_acc, org_name]]
		if org_name in name2final_ver:
			name2final_ver[org_name].append(asm_acc)
		else:
			name2final_ver[org_name] = [asm_acc]

	# now keep only the final version: GCF > GCA if avail, most recent ver
	for name in name2final_ver:
		name2final_ver[name] = sorted(name2final_ver[name])[-1]
	return taxid2asmnames, name2final_ver


# Given organism names for each taxid, assign each a unique taxid
def assign_unique_taxids(taxid2asmnames, name2final_ver):
	asm2uniq_taxid = {}
	for taxid in taxid2asmnames:
		# query name2final_ver for each assembly with this taxid,
		#  	only keep final versions. then assign unique taxids for each kept.
		final_names = [i[0] for i in taxid2asmnames[taxid]
							if i[0] == name2final_ver[i[1]]]
		if len(final_names) == 1:
			asm2uniq_taxid[final_names[0]] = taxid
		else:
			for i in range(len(final_names)):
				asm2uniq_taxid[final_names[i]] = taxid + '.' + str(i)
	return asm2uniq_taxid


# Given assembly names to taxids, use taxtree and genomic file to extract rest
#  	of needed info (organism length and lineage), & uncompress genomic file
def build_dbinfo_and_extract(args, asm2taxid, taxtree):
	dbinfo = open(args.output_dir + 'db_info.txt', 'w')
	dbinfo.write('Accesion\tLength   \tTaxID   \tLineage \tTaxID_Lineage\n')
	dbinfo.write('Unmapped\t0\tUnmapped\t|||||||Unmapped\t|||||||Unmapped\n')
	for asm_acc in asm2taxid:
		# for each assembly, generate its lineage info and output file handler
		taxid, genome_len, acc_list = asm2taxid[asm_acc], 0, []
		name_lin, taxid_lin = trace_lineages(taxid.split('.')[0], taxtree)
		if '.' in taxid:  # need to place unique taxid into taxidlin
			taxid_lin = '|'.join(taxid_lin.split('|')[:-1] + [taxid])
		outname = 'taxid_' + taxid.replace('.', '_') + '_genomic.fna'

		# parse the genomic file for genome length and contained accessions
		with(gzip.open(args.input_dir+asm_acc+'_genomic.fna.gz','rt')) as infile:
			with(open(args.output_dir + outname, 'w')) as outfile:
				for line in infile:
					outfile.write(line)  # write uncompressed line to output loc
					if line.startswith('>'):
						acc_list.append(line.split()[0][1:])
					else:
						genome_len += len(line.strip())

		# write the gathered information to the db_info file
		for acc in acc_list:
			dbinfo.write('\t'.join([acc, str(genome_len), taxid,
										name_lin, taxid_lin]) + '\n')
	dbinfo.close()


def main():
	args = parseargs()
	if not args.input_dir.endswith('/'):
		args.input_dir += '/'
	if not args.output_dir.endswith('/'):
		args.output_dir += '/'
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	# get taxonomic information
	taxtree, taxid2namelin = get_taxonomy_info()

	# gather microbial assembly names, then assign them to unique taxids
	taxid2asmnames, name2final_ver = get_taxids2asmnames(args, taxid2namelin)
	asm2uniq_taxid = assign_unique_taxids(taxid2asmnames, name2final_ver)
	# now build dbinfo and gunzip genome files to args.output_dir
	build_dbinfo_and_extract(args, asm2uniq_taxid, taxtree)


if __name__ == '__main__':
	main()
#
