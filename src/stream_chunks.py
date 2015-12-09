#!/usr/bin/python

import Bio.SeqIO, argparse, gzip, math, sys

def parse_arguments():
	""" Parse commandline arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="""Stream sequence from multi-FASTA file in defined chunks""")
	parser.add_argument('--fasta', type=str, required=True,
						help="Path to input multi-FASTA file")
	parser.add_argument('--chunk_size', type=int, required=False,
						help="Chunk size in base-pairs")
	args = vars(parser.parse_args())
	return args

def stream_chunks(args):
	ext = args['fasta'].split('.')[-1]
	infile = gzip.open(args['fasta']) if ext == 'gz' else open(args['fasta'])
	for rec in Bio.SeqIO.parse(infile, 'fasta'):

		chunks = int(math.ceil(len(rec.seq)/float(args['chunk_size'])))
		for chunk in range(chunks):
			start = (chunk + 1) * args['chunk_size'] - args['chunk_size']
			end = (chunk + 1) * args['chunk_size']
			seq = rec.seq[start:end].upper()
			sys.stdout.write('>%s_%s\n%s\n' % (rec.id, chunk, seq))

if __name__ == "__main__":
	args = parse_arguments()
	stream_chunks(args)