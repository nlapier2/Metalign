# takes a directory path as input, combines CAMI profiles in that directory, and
#  ensures compatibility with fingerprinting
import glob, sys
indir, outname = sys.argv[1], sys.argv[2]
if not indir.endswith('/'):
	indir += '/'

with(open(outname, 'w')) as outfile:
	for fname in glob.glob(indir + '*'):
		with(open(fname, 'r')) as infile:
			for line in infile:
				if line.startswith('@SampleID') or line.startswith('@Ranks') or len(line) < 5:
					outfile.write(line)
				elif line.startswith('@Version'):
					outfile.write('@Version:0.9\n')
				elif line.startswith('@@TAXID'):
					outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n')
				else:
					splits = line.strip().split('\t')
					if splits[0].count('.') == 2:
						splits[0] = splits[0].rsplit('.', 1)[0] + '0' + splits[0].rsplit('.', 1)[1]
					if splits[2].count('.') == 2:
						splits[2] = splits[2].rsplit('.', 1)[0] + '0' + splits[2].rsplit('.', 1)[1]
					outfile.write('\t'.join(splits[:5]) + '\n')
#
