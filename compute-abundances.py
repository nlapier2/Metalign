import argparse, operator, os, random, sys, time


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
virus_info = __location__ + '/data/accession2info-viral.txt'
fungi_info = __location__ + '/data/accession2info-fungi.txt'
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
arg_defaults = {'abundance_cutoff': -1, 'min_map': -1, 'max_ed': -1, 'pct_id': -1, 'read_cutoff': -1}
fun_defaults = {'abundance_cutoff': -1, 'min_map': 100, 'max_ed': 1, 'pct_id': 0.99, 'read_cutoff': 100}
vir_defaults = {'abundance_cutoff': -1, 'min_map': 0, 'max_ed': 999999, 'pct_id': 0.6, 'read_cutoff': 10}


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Compute abundance estimations for species in a sample.')
    parser.add_argument('sam', help='.sam file or file with list of SAM files to process. Required.')
    parser.add_argument('--abundance_cutoff', type=float, default=-1, help='Organism abundance to count it as present.')
    parser.add_argument('--fungi', action='store_true', help='Profile fungi.')
    parser.add_argument('--min_map', type=int, default=-1, help='Minimum bases mapped to count a hit.')
    parser.add_argument('--max_ed', type=int, default=-1, help='Maximum edit distance from a reference to count a hit.')
    parser.add_argument('--normalize', type=bool, choices=[True,False], default=True,
    					help='Normalize species abundance by genome length or not. Default: True')
    parser.add_argument('--output', default='abundances.txt', help='Output abundances file. Default: abundances.txt')
    #parser.add_argument('--paired', action='store_true', default=False, help='Use if reads are paired end.')
    parser.add_argument('--pct_id', type=float, default=-1, help='Minimum percent identity from reference to count a hit.')
    parser.add_argument('--read_cutoff', type=int, default=-1, help='Number of reads to count an organism as present.')
    parser.add_argument('--virus', action='store_true', help='Profile viruses.')
    args = parser.parse_args()
    return args


def set_params(args):
	if args.virus:
		if args.abundance_cutoff == arg_defaults['abundance_cutoff']:
			args.abundance_cutoff = vir_defaults['abundance_cutoff']
		if args.min_map == arg_defaults['min_map']:
			args.min_map = vir_defaults['min_map']
		if args.max_ed == arg_defaults['max_ed']:
			args.max_ed = vir_defaults['max_ed']
		if args.pct_id == arg_defaults['pct_id']:
			args.pct_id = vir_defaults['pct_id']
		if args.read_cutoff == arg_defaults['read_cutoff']:
			args.read_cutoff = vir_defaults['read_cutoff']
	else:
		if args.abundance_cutoff == arg_defaults['abundance_cutoff']:
			args.abundance_cutoff = fun_defaults['abundance_cutoff']
		if args.min_map == arg_defaults['min_map']:
			args.min_map = fun_defaults['min_map']
		if args.max_ed == arg_defaults['max_ed']:
			args.max_ed = fun_defaults['max_ed']
		if args.pct_id == arg_defaults['pct_id']:
			args.pct_id = fun_defaults['pct_id']
		if args.read_cutoff == arg_defaults['read_cutoff']:
			args.read_cutoff = fun_defaults['read_cutoff']
	return args


def find_taxid(tag):
        if not '|' in tag:
                return tag
        else:
                splits = tag.split('|')
                for sp in splits:
                        if '_' in sp:
                                return sp
        return tag


def ids2info(args):
	acc2info, clade2gi, lin2len = {}, {}, {}
	if args.virus:
		infofile = virus_info
	else:
		infofile = fungi_info
	with(open(infofile, 'r')) as infile:
		for line in infile:
			splits = line.strip().split('\t')
			if len(splits) == 4:
				acc2info[splits[0]] = [splits[1], splits[2], splits[3]]
				if splits[3] in lin2len:
					lin2len[splits[3]] += float(splits[1])
				else:
					lin2len[splits[3]] = float(splits[1])
			else:
				clade2gi[splits[1]] = splits[0]
	return acc2info, clade2gi, lin2len


def filter_line(args, splits):  # determine whether to filter this line out
    cigar = splits[5]
    if cigar == '*':
    	return True
    matched_len, total_len, cur = 0, 0, 0
    for ch in cigar:
    	if not ch.isalpha():  # ch is a number
    		cur = (cur * 10) + int(ch)
    	else:  # ch is a letter
    		if ch == 'M' or ch == '=':
    			matched_len += cur
    		total_len += cur
    		cur = 0
    if matched_len < args.min_map:
    	return True
    edit_distance = int(splits[11][5:])
    if edit_distance > args.max_ed:
    	return True
    #if float(mapped) / float(mapped + edit_distance) < args.pct_id:
    if float(matched_len) / float(total_len) < args.pct_id:
    	return True
    return False  # if read passes quality checks, don't filter


def compute_abundances(args, samfile, acc2info, clade2gi, lin2len):
	infile = open(samfile, 'r')
	ids, ids2abs, clade2abs = [], {}, {}
	prev_read_num, prev_tag, prev_count, ignore = '', '', 1.0, False
	multimapped, ids2reads, read_ordering = {}, {}, []
	lc = 0

	print 'Reading sam file ' + samfile
	for line in infile:
		lc += 1
		if lc % 1000000 == 0:
			print 'Done reading ' + str(lc) + ' lines of sam file'
		if line.startswith('@'):
			continue
		splits = line.strip().split('\t')
		if filter_line(args, splits):
			continue

		tag = find_taxid(splits[2])
		if tag == '*':
			continue
		#read_num = int(splits[0])
		read_num = splits[0]
		read_ordering.append(read_num)
		if read_num == prev_read_num and tag == prev_tag:
			pass
			#if prev_count < 2.0 and args.paired == True:
			#	prev_count += 1.0
		elif read_num == prev_read_num and tag != prev_tag:
			ignore = True
			strnum = str(prev_read_num)
			if strnum not in multimapped:
				multimapped[strnum] = [prev_tag]
			else:
				multimapped[strnum].append(prev_tag)
			prev_tag = tag
		else:
			if not(prev_read_num == '' or ignore == True):
				if prev_tag not in ids2reads:
					ids2reads[prev_tag] = [prev_read_num]
				else:
					ids2reads[prev_tag].append(prev_read_num)
				if prev_tag not in ids:
					ids.append(prev_tag)
				if prev_tag in ids2abs:
					ids2abs[prev_tag] += prev_count
				else:
					ids2abs[prev_tag] = prev_count
			elif ignore == True:
				multimapped[prev_read_num].append(prev_tag)
	                prev_read_num = read_num
	                prev_tag = tag
	                prev_count = 1.0
	                ignore = False
	infile.close()
	print 'Done reading sam file.'

	if not(prev_read_num == '' or ignore == True):
		if prev_tag not in ids2reads:
	                ids2reads[prev_tag] = [prev_read_num]
	        else:
	                ids2reads[prev_tag].append(prev_read_num)
	        if prev_tag not in ids:
	                ids.append(prev_tag)
	       	if prev_tag in ids2abs:
	                ids2abs[prev_tag] += prev_count
	       	else:
	                ids2abs[prev_tag] = prev_count
	elif ignore == True:
		multimapped[read_num].append(prev_tag)

	print 'Deleting clades with insufficient evidence...'
	clade2ids, del_list = {}, []
	for taxid in ids2abs.keys():
		clade = acc2info[taxid][2]
		if clade in clade2abs:
			clade2abs[clade] += ids2abs[taxid]
			clade2ids[clade].append(taxid)
		else:
			clade2abs[clade] = ids2abs[taxid]
			clade2ids[clade] = [taxid]
	for clade in clade2abs.keys():
		if clade2abs[clade] < args.read_cutoff:
			del_list.append(clade)  # mark for deletion due to insufficient evidence
		else:
			clade2abs[clade] /= lin2len[clade]  # normalize
	for key in del_list:
		for taxid in clade2ids[key]:
			del ids2abs[taxid]
		del clade2abs[key]
	print 'Done deleting clades.'

	print 'Assigning multimapped reads...'
	added = {}
	for read in multimapped.keys():
		randnum = random.random()
		options, total = {}, 0.0
		for taxid in multimapped[read]:
			if taxid not in ids2abs or taxid in options:
				continue
			ab = ids2abs[taxid]
			total += ab
			options[taxid] = ab
		if len(options.keys()) == 0:
			continue
		for key in options.keys():
			ab = options[key] / total
			if ab >= randnum:
				if key in added:
					added[key] += 1.0 #2.0
				else:
					added[key] = 1.0 #2.0
				break
			else:
				randnum -= ab

	for key in added.keys():
		clade = acc2info[key][2]
		clade2abs[clade] += (added[key] / lin2len[clade])
	print 'Multimapped reads assigned.'

	total_ab = 0.0
	for clade in clade2abs.keys():
		total_ab += clade2abs[clade]

	for clade in clade2abs.keys():
		clade2abs[clade] = clade2abs[clade] * 100.0 / total_ab  # normalize abundances

	results = {}  # holds results for cami format
	for clade in clade2abs:
		all_taxa = clade.split('|')
		all_gi = [clade2gi[i] if i in clade2gi else '?' for i in all_taxa]
		results[clade] = [all_gi[-1], 'strain', '|'.join(all_gi), clade, clade2abs[clade]]

		all_levels = ['|'.join(all_taxa[:i]) for i in range(1, 1+len(all_taxa))]
		for i in range(len(all_levels)-1):
			if all_levels[i] in results:
				results[all_levels[i]][-1] += clade2abs[clade]
			else:
				results[all_levels[i]] = [all_gi[i], RANKS[i],'|'.join(all_gi[:i+1]),all_levels[i],clade2abs[clade]]

	return results


def main():
	args = set_params(parseargs())
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		print 'Error: --pct_id must be between 0.0 and 1.0, inclusive'
		sys.exit()
	if (args.virus and args.fungi) or not(args.virus or args.fungi):
		print 'Error: must specify either --virus or --fungi.'
		sys.exit()
	samfiles = []
	if args.sam.endswith('.sam'):
		samfiles.append(args.sam)
	else:
		with(open(args.sam, 'r')) as filenames:
			for line in filenames:
				samfiles.append(line.strip())
	acc2info, clade2gi, lin2len = ids2info(args)  # maps NCBI taxID to length, GI, lineage

	results = {}
	for sam in samfiles:
		res = compute_abundances(args, sam, acc2info, clade2gi, lin2len)
		for clade in res:
			if clade not in results:
				results[clade] = res[clade]
			else:
				results[clade][-1] += res[clade][-1]

	lev_res = {i:[] for i in range(len(RANKS))}
	for clade in results:
		results[clade][-1] /= len(samfiles)  # average over all input files
		lev_res[len(clade.split('|'))-1].append(results[clade])  # accumulate results by tax level

	print 'Writing clade abundances...'
	with(open(args.output, 'w')) as outfile:
		outfile.write('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n')
		for i in range(len(RANKS)):
			lines = lev_res[i]
			lines.sort(key=lambda x: 100.0-x[-1])
			if lines == None or len(lines) < 1:
				continue
			for line in lines:
				line = [str(i) for i in line]
				outfile.write('\t'.join(line)+'\n')
	print 'Done.'


if __name__ == '__main__':
	main()
#
