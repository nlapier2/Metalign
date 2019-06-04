import argparse, glob, gzip, os, subprocess


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given refseq" +
		"database and taxonomy files, build Metalign db.")
	parser.add_argument('--input_file', required = True,
		help = 'Name of refseq database file.')
	parser.add_argument('--map_file', required = True,
		help = 'Name of file mapping accession to taxID.')
	parser.add_argument('--taxonomy_dir', required = True,
		help = 'Name of directory with names.dmp and nodes.dmp.')
	parser.add_argument('--output_dir', default='organism_files/',
		help = 'Name of directory to write selected genome files to.')
	parser.add_argument('--no_output', action = 'store_true',
		help = 'Do not split; just (re-)write db_info.')
	args = parser.parse_args()
	return args


def build_taxtree(taxonomy_dir):
	# taxtree is stored as a dict with TaxID keys linking to
	# 		taxonomic name, taxonomic rank, and parent TaxID
	# Thus, it will be easy to traverse through lineage given starting TaxID
	taxtree = {}

	with(open(taxonomy_dir + 'names.dmp', 'r')) as names:
		for line in names:
			if 'scientific name' not in line:
				continue  # skip alternate names and other unnecessary info
			taxid = line.split()[0]
			name = line.split('|')[1].strip()
			taxtree[taxid] = [name]  # Key is TaxID, value is tax. name

	with(open(taxonomy_dir + 'nodes.dmp', 'r')) as nodes:
		for line in nodes:
			# add taxonomic rank and parent TaxID to list of info for this taxID
			splits = [field.strip() for field in line.split('|')]
			taxid, parent, rank = splits[:3]
			taxtree[taxid].extend([rank, parent])
			#taxtree[splits[0]].extend([splits[4], splits[2]])
	return taxtree


# Uses taxtree to find lineage of taxid in both taxonomic names and IDs
def trace_lineages(taxid, taxtree):
	cami_ranks = {'superkingdom': 0, 'phylum': 1, 'class': 2, 'order': 3,
		'family': 4, 'genus': 5, 'species': 6, 'strain': 7}
	name_lineage, taxid_lineage = ['' for i in range(8)], ['' for i in range(8)]
	cur_taxid = taxid

	# Record lowest level taxonomic info if it's below species level
	if cur_taxid in taxtree:
		name, rank, parent = taxtree[cur_taxid]
	else:
		return 'NONE', 'NONE'
	if rank not in cami_ranks:  # strains are recorded as 'no rank' in nodes.dmp
		name_lineage[-1] = name
		taxid_lineage[-1] = cur_taxid
		cur_taxid = parent  # traverse up to parent taxid

	# Now traverse up the tree
	while cur_taxid != '1':  # while we're not at the root
		if cur_taxid in taxtree:
			name, rank, parent = taxtree[cur_taxid]
		else:
			return 'NONE', 'NONE'
		if rank in cami_ranks:  # if this is a relevant CAMI taxonomic rank
			index = cami_ranks[rank]  # then record this node's info at the rank
			name_lineage[index] = name
			taxid_lineage[index] = cur_taxid
		cur_taxid = parent  # traverse up to parent taxid and repeat
	return '|'.join(name_lineage), '|'.join(taxid_lineage)


def parse_map(map_file):
	acc2taxid = {}
	with(open(map_file, 'r')) as infile:
		for line in infile:
			gi_num, acc, taxid = line.strip().split()
			acc2taxid[acc] = taxid
	return acc2taxid


def split_and_process(input_file, output_dir, taxtree, acc2taxid, block_output):
	outfile, acc2info, taxids_processed = '', {}, {}
	suspend = False
	unwanted = ['Chordata', 'Streptophyta', 'Arthropoda', 'Echinodermata', 'Mollusca']
	with(open(input_file, 'r')) as infile:
		acc, final_acc, taxid, name_lin, tax_lin = '', '', '', '', ''
		for line in infile:
			if line.startswith('>'):
				# get the new accession and its taxid and lineage information
				acc, name = line.strip().split('|')[3:]
				acc_prefix = acc[:7]  # heuristic for same strain
				name = name.strip().split(',')[0]
				taxid = acc2taxid[acc]
				name_lin, tax_lin = trace_lineages(taxid, taxtree)

				# if this is an unwanted taxon, suspend output
				suspend = False
				for unwanted_taxon in unwanted:
					if unwanted_taxon in name_lin:
						suspend = True
						break
				if 'plasmid' in name or 'Plasmid' in name:
					suspend = True
				if suspend:
					continue

				# prevent duplicate accesions/taxids by adding '.X' with unique num. X
				final_acc, version = acc, 0
				while final_acc in acc2info:
					version += 1
					final_acc = acc + '.' + str(version)
				final_tax, version = taxid, 0
				while final_tax in taxids_processed:
					if acc_prefix == taxids_processed[final_tax]:
						break
					version += 1
					final_tax = taxid + '.' + str(version)
				if final_tax != taxid:
					namesplits, taxsplits = name_lin.split('|'), tax_lin.split('|')
					taxsplits[-1] = final_tax
					namesplits[-1] = name
					name_lin, tax_lin = '|'.join(namesplits), '|'.join(taxsplits)

				line = line.replace(acc, final_acc)
				acc2info[final_acc] = [0, final_tax, name_lin, tax_lin]
				taxids_processed[final_tax] = acc_prefix

				if not block_output:
					#outfile = open(output_dir + 'acc_' + acc.replace('.', '-') + '_genomic.fna', 'w')
					outfile = open(output_dir + 'taxid_' + final_tax.replace(
						'.', '_') + '_genomic.fna', 'a')
			elif not suspend:
				acc2info[final_acc][0] += len(line.strip())
			if not block_output and not suspend:
				outfile.write(line)  # write this line to the taxid genomic file
	return acc2info


def write_dbinfo(output_dir, acc2info):
	with(open(output_dir + 'db_info.txt', 'w')) as outfile:
		outfile.write('Accesion\tLength   \tTaxID   \tLineage \tTaxID_Lineage\n')
		outfile.write('Unmapped\t0\tUnmapped\t|||||||Unmapped\t|||||||Unmapped\n')
		for acc in acc2info:
			acc2info[acc][0] = str(acc2info[acc][0])
			outfile.write(acc + '\t' + '\t'.join(acc2info[acc]) + '\n')


def main():
	args = parseargs()
	if not args.taxonomy_dir.endswith('/'):
		args.taxonomy_dir += '/'
	if not args.output_dir.endswith('/'):
		args.output_dir += '/'
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)

	taxtree = build_taxtree(args.taxonomy_dir)
	acc2taxid = parse_map(args.map_file)
	acc2info = split_and_process(args.input_file, args.output_dir, taxtree, acc2taxid, args.no_output)
	write_dbinfo(args.output_dir, acc2info)


if __name__ == '__main__':
	main()
#
