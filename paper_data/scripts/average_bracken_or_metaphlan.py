import argparse, glob


def parseargs():
	parser = argparse.ArgumentParser(description="Given filename pattern" +
				" for bracken or metaphlan input files.")
	parser.add_argument('--pattern', required=True,
		help='Pattern to match, e.g. "/path/to/profiles/*pattern*". Required.')
	parser.add_argument('--output', required=True,
		help='Where to write output profile. Required.')
	args = parser.parse_args()
	return args


def parse_profile(profile):
	profile_results = {}
	with(open(profile, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.strip().split('\t')
			taxon, abundance = splits[0], float(splits[-1])
			profile_results[taxon] = abundance
	return profile_results


# Writes results out in simple format
def write_results(output, results):
	with(open(output, 'w')) as outfile:
		for taxon in results:
			outfile.write(taxon + '\t' + str(results[taxon]) + '\n')


def main():
	args = parseargs()
	profiles = glob.glob(args.pattern)
	num_profiles = float(len(profiles))
	results = {}

	# get results for each profile and sum them up
	for pro in profiles:
		profile_results = parse_profile(pro)
		for taxon in profile_results:
			averaged_abundance = profile_results[taxon] / num_profiles
			if taxon not in results:
				results[taxon] = averaged_abundance
			else:
				results[taxon] += averaged_abundance
	write_results(args.output, results)


if __name__ == '__main__':
	main()
#
