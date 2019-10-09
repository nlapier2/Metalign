# Takes in a CAMI format file and returns a text file formatted for Krona
import argparse


def parseargs():
	parser = argparse.ArgumentParser(description="Given CAMI format input " +
				"file, generate a text file that Krona can make a chart from.")
	parser.add_argument('--input', required=True,
		help='CAMI-format input file. Required.')
	parser.add_argument('--output', required=True,
		help='Where to write output text file. Required.')
	parser.add_argument('--no_strains', action='store_true',
		help='Use if strains are not recorded for this file.')
	args = parser.parse_args()
	return args


def main():
	args = parseargs()
	with(open(args.input, 'r')) as infile:
		with(open(args.output, 'w')) as outfile:
			for line in infile:
				if line.startswith('@') or line.startswith('#') or len(line)<5:
					continue  # skip headers and blank lines
				if args.no_strains and 'species' not in line:
					continue
				if not args.no_strains and 'strain' not in line:  # for Krona, only lowest level needed
					continue
				splits = line.strip().split('\t')
				namelin, abundance = splits[3].split('|'), splits[4]
				namelin = [i if i != '' else 'unlabeled taxon' for i in namelin]
				outfile.write(abundance + '\t' + '\t'.join(namelin) + '\n')


if __name__ == '__main__':
	main()
#
