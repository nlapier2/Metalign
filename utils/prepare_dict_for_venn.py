# Takes in CAMI and/or MetaPhlAn2 format file(s) as well as a taxonomic level,
#    and writes a dict with method mapped to taxa at that level, used to make venn diagram
import argparse, sys


def parseargs():
	parser = argparse.ArgumentParser(description="Given CAMI/MetaPhlAn2 format input files " +
				"as well as a taxonomic level, prepare dict used to make venn diagram.")
	parser.add_argument('--output', required=True,
		help='Where to write output text file. Required.')
	parser.add_argument('--tax_level', required=True,
		choices = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
		help='Which taxonomic level to extract.')
	parser.add_argument('--cami_input', nargs='+',
		help='CAMI-format input file(s).')
	parser.add_argument('--metaphlan_input', nargs='+',
		help='MetaPhlAn2-format input file(s).')
	parser.add_argument('--cami_names', nargs='+',
		help='Method names corresponding to CAMI input files.')
	parser.add_argument('--metaphlan_names', nargs='+',
		help='Method names corresponding to MetaPhlAn2 input files.')
	parser.add_argument('--abundance_cutoff', type=float, default=0.01,
		help='Do not include taxa with below this abundance threshold.')
	args = parser.parse_args()
	return args


def get_set_from_cami_file(cami_file, tax_level, abundance_cutoff):
	taxa_list = []
	with(open(cami_file, 'r')) as infile:
		for line in infile:
			if line.startswith('#') or line.startswith('@') or len(line) < 5:
				continue
			splits = line.split('\t')
			if splits[1] != tax_level:
				continue
			if float(splits[4]) < abundance_cutoff:
				continue
			namelin = splits[3]
			name = namelin.split('|')[-1]
			taxa_list.append(name)
	taxa_list.sort()
	return set(taxa_list)


def translate_taxlevel_to_mp(tax_level):
	if tax_level == 'superkingdom': return 'k'
	if tax_level == 'phylum': return 'p'
	if tax_level == 'class': return 'c'
	if tax_level == 'order': return 'o'
	if tax_level == 'family': return 'f'
	if tax_level == 'genus': return 'g'
	if tax_level == 'species': return 's'
	if tax_level == 'strain': return 't'


def get_set_from_metaphlan_file(metaphlan_file, tax_level, abundance_cutoff):
	tax_level = translate_taxlevel_to_mp(tax_level)
	taxa_list = []
	with(open(metaphlan_file, 'r')) as infile:
		for line in infile:
			if len(line) < 5:
				continue
			splits = line.strip().split('\t')
			taxlin, abundance = splits[0], float(splits[1])
			if abundance < abundance_cutoff:
				continue
			taxon = taxlin.split('|')[-1]
			if '__' not in taxon:  # this is bracken's almost metaphlan format
				taxon = taxon.replace('_', '__')
				if taxon[0] == 'd':
					taxon = 'k' + taxon[1:]
			splits = taxon.split('__')
			this_taxlevel, this_name = splits[0], splits[1]
			if this_taxlevel != tax_level:
				continue
			taxa_list.append(this_name)
	taxa_list.sort()
	return set(taxa_list)


def main():
	args = parseargs()
	if len(args.cami_names) != len(args.cami_input):
		sys.exit("Error: number of --cami_names should equal total number of --cami_input files.")
	if len(args.metaphlan_names) != len(args.metaphlan_input):
		sys.exit("Error: number of --metaphlan_names should equal number of --metaphlan_input files.")

	method_to_taxa = {}
	for i in range(len(args.cami_input)):
		taxa_set = get_set_from_cami_file(
			args.cami_input[i], args.tax_level, args.abundance_cutoff)
		method_to_taxa[args.cami_names[i]] = taxa_set
	for i in range(len(args.metaphlan_input)):
		taxa_set = get_set_from_metaphlan_file(
			args.metaphlan_input[i], args.tax_level, args.abundance_cutoff)
		method_to_taxa[args.metaphlan_names[i]] = taxa_set

	with(open(args.output, 'w')) as outfile:
		outfile.write(str(method_to_taxa))


if __name__ == '__main__':
	main()
#
