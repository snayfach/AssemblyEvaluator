#!/usr/bin/python

import os, argparse, gzip, sys, subprocess

def parse_arguments():
	""" Parse commandline arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="""Evaluate a metagenomic assembly by mapping to a reference genome""")
	parser.add_argument('--query', type=str, required=True,
						help="Path to input multi-FASTA file containing one or more contigs")
	parser.add_argument('--ref', type=str, required=False,
						help="Path to reference genome FASTA file")
	parser.add_argument('--out', dest='out', type=str, required=True,
						help="""
						Output directory:
						<out>/alignments.m8 (tabular alignments)
						""")
	parser.add_argument('--threads', type=int, default=1,
						help="Threads to use for BLAST")
	parser.add_argument('--chunk_size', type=int, default=1000,
						help="Chunk size (in bp) for splitting up queries")
	args = vars(parser.parse_args())
	if not os.path.isdir(args['out']): os.mkdir(args['out'])
	add_binaries(args)
	return args

def add_binaries(args):
	""" Add paths to binaries """
	script_path = os.path.abspath(__file__)
	script_dir = os.path.dirname(script_path)
	main_dir = os.path.dirname(script_dir)
	args['blastp'] = '%s/bin/blastp' % main_dir
	args['blastn'] = '%s/bin/blastn' % main_dir
	args['makeblastdb'] = '%s/bin/makeblastdb' % main_dir
	args['stream_chunks'] = '%s/stream_chunks.py' % script_dir

def run_shell_command(command):
	""" Capture stdout, stderr. Check unix exit code and exit if non-zero """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s" % (command, err)
		sys.exit(err_message)
	else:
		return out

def run_blast(args):
	""" Run blastn """
	cmd = '%s ' % args['stream_chunks']
	cmd += '--fasta %s ' % args['query']
	cmd += '--chunk_size %s | ' % args['chunk_size']
	cmd += '%s ' % args['blastn']
	cmd += '-query /dev/stdin '
	cmd += '-subject %s ' % args['ref']
	cmd += '-outfmt 6 '
	cmd += '-max_target_seqs 1 '
	cmd += '-num_threads %s ' % args['threads']
	out = run_shell_command(cmd)
	return out

def parse_blast(m8_path):
	""" Yield formatted record from BLAST m8 file """
	formats = {
		0:str, 1:str, 2:float, 3:int, 4:float,
		5:float, 6:float, 7:float, 8:float,
		9:float, 10:float, 11:float
	}
	fields = {
		0:'query_id',1:'target_id',2:'pid',3:'aln',4:'mis',
		5:'gaps',6:'qstart',7:'qend',8:'tstart',
		9:'tend',10:'evalue',11:'score'
	}
	for line in open(m8_path):
		r = line.rstrip().split()
		if len(r) != 12: continue
		else: yield dict([(fields[i], formats[i](v)) for i, v in enumerate(r)])

def add_aln_cov(m8, queries):
	""" Add query and target coverage to m8 record """
	qlen = float(queries[m8['query_id']].length)
	tlen = float(m8['target_id'].split('_')[-3])
	m8['qcov'] = round(m8['aln']/qlen, 2)
	m8['tcov'] = round(m8['aln']/tlen, 2)
	return m8

def write_m8(args, blastout):
	""" Write tab delimited m8 file to disk """
	outfile = open(os.path.join(args['out'], 'alignments.m8'), 'w')
	outfile.write(blastout)
	outfile.close()


if __name__ == "__main__":

	args = parse_arguments()
	blastout = run_blast(args)
	write_m8(args, blastout)




