import argparse, os, shlex, subprocess, sys

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
bwa_exec = __location__ + '/bwa.kit/bwa'
virus_ind = __location__ + '/data/viral-refseq.fna'
fungi_ind = __location__ + '/data/fungi-refseq.fna'


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Run BWA mem for viruses or fungi.')
    parser.add_argument('reads', nargs='+', help='Reads to run BWA on. One or two files, depending on single or paired end. Required.')
    parser.add_argument('--fungi', action='store_true', help='Run BWA on fungi. This or --virus required.')
    parser.add_argument('--virus', action='store_true', help='Run BWA on viruses. This or --fungi required.')
    parser.add_argument('--output', default='alignments.sam', help='Where to output BWA results to.')
    parser.add_argument('--paired', action='store_true', help='Treat input reads file as interleaved paired end reads.')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if (args.virus and args.fungi) or not(args.virus or args.fungi):
        print 'Error: must specify either --virus or --fungi.'
        sys.exit()
    if len(args.reads) > 2:
        print 'Error: must specify one (single end) or two (paired end) reads files.'
        sys.exit()

    args.reads = ' '.join(args.reads)
    paired = '-p' if len(args.reads) == 1 and args.paired else ''

    if args.virus:
        cmd = ' '.join([bwa_exec, 'mem', paired, '-a', virus_ind, args.reads])
    else:  # args.fungi
        cmd = ' '.join([bwa_exec, 'mem', paired, '-a', fungi_ind, args.reads])
    with(open(args.output, 'w')) as outfile:
        subprocess.call(shlex.split(cmd), stdout=outfile)


if __name__ == '__main__':
	main()
#
