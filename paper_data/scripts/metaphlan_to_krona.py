# Takes in a metaphlan2 file FROM MY AVERAGING SCRIPT and puts in format for Krona
import argparse


def parseargs():
	parser = argparse.ArgumentParser(description="Generate a text file " +
				"that Krona can make a chart from.")
	parser.add_argument('--input', required=True,
		help='Input file. Required.')
	parser.add_argument('--output', required=True,
		help='Where to write output text file. Required.')
	parser.add_argument('--from_bracken', action='store_true',
		help='Use if --input was generated by Bracken helper script.')
	args = parser.parse_args()
	return args


def main():
	args = parseargs()
	with(open(args.input, 'r')) as infile:
		with(open(args.output, 'w')) as outfile:
			for line in infile:
				if not args.from_bracken and 't__' not in line:
					continue
				if args.from_bracken and 's_' not in line:
					continue
				splits = line.strip().split('\t')
				namelin, abundance = splits[0], splits[1]
				if not args.from_bracken:
					namelin = namelin.replace('__', '').split('|')
				else:
					namelin = namelin.replace('_', '').split('|')
				namelin = [i[1:] for i in namelin]
				outfile.write(abundance + '\t' + '\t'.join(namelin) + '\n')


if __name__ == '__main__':
	main()
#