# Takes in a filename pattern matching two or more CAMI-format metagenome
#  	profiles, averages the taxa abundances across all files, and writes results
import argparse, glob


RANKS = ['superkingdom', 'phylum', 'class', 'order',
		'family', 'genus', 'species', 'strain']


def parseargs():
	parser = argparse.ArgumentParser(description="Given filename pattern" +
				" for CAMI-format profiles, average abundances & output.")
	parser.add_argument('--pattern', required=True,
		help='Pattern to match, e.g. "/path/to/profiles/*pattern*". Required.')
	parser.add_argument('--output', required=True,
		help='Where to write output profile. Required.')
	parser.add_argument('--sampleID', default='NONE',
		help='Sample ID for output. Defaults to the value for --pattern.')
	args = parser.parse_args()
	return args


def parse_profile(profile):
	profile_results = {}
	with(open(profile, 'r')) as infile:
		for line in infile:
			if line.startswith('@') or line.startswith('#') or len(line) < 5:
				continue  # skip headers and blank lines
			splits = line.strip().split('\t')
			taxid = splits[0]
			splits[4] = float(splits[4])  # abundance, which we want as a float
			profile_results[taxid] = splits
	return profile_results


# groups results into the different taxonomic ranks for outputting purposes
def generate_rankwise_results(results, num_profiles):
	rank_results = {}
	for i in range(len(RANKS)):
		rank_results[i] = []
	for taxon in results:
		results[taxon][4] = results[taxon][4] / num_profiles  #avg over profiles
		rank = RANKS.index(results[taxon][1])
		rank_results[rank].append(results[taxon])  # results per rank
	return rank_results


# Writes results out in CAMI format
def write_results(output, sampleID, rank_results):
	with(open(output, 'w')) as outfile:
		# Print some CAMI format header lines
		outfile.write('@SampleID:' + sampleID + '\n')
		outfile.write('@Version:Metalign-v0.2\n')
		outfile.write('@Ranks: ' +
			'superkingdom|phylum|class|order|family|genus|species|strain\n\n')
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\t' +
			'PERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n')

		for i in range(len(RANKS)):
			lines = rank_results[i]  # all lines to write for this tax level
			# now sort clades in rank by descending abundance
			lines.sort(key=lambda x: 100.0-x[4])
			if lines == None:
				continue
			for line in lines:
				if line[4] < 0.00001:
					line[4] = 0.00001
				else:
					line[4] = float('%.5f' % line[4])
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')


def main():
	args = parseargs()
	if args.sampleID == 'NONE':
		args.sampleID = args.pattern  # sampleID defaults to pattern name
	profiles = glob.glob(args.pattern)
	num_profiles = float(len(profiles))
	results = {}

	# get results for each profile and sum them up
	for pro in profiles:
		profile_results = parse_profile(pro)
		for taxon in profile_results:
			if taxon not in results:
				results[taxon] = profile_results[taxon]
			else:
				results[taxon][4] += profile_results[taxon][4]

	rankwise_results = generate_rankwise_results(results, num_profiles)
	write_results(args.output, args.sampleID, rankwise_results)


if __name__ == '__main__':
	main()
#
