#!/usr/bin/python

import os, argparse, gzip, sys, subprocess, Bio.SeqIO

##
## Classes
##

class Genome:
	""" Base class for genome sequences """
	def __init__(self, args):
		self.name = os.path.basename(args['genome'])
		self.contigs = self.init_contigs(args)
		self.count_contigs = len(self.contigs.keys())
		self.length = self.length()
		self.uncalled = self.uncalled_bp()
		self.length_trim = self.length - self.uncalled
		
	def init_contigs(self, args):
		contigs = {}
		for rec in Bio.SeqIO.parse(args['genome'], 'fasta'):
			contigs[rec.id] = Contig(rec)
		return contigs

	def length(self):
		return(sum([c.length for c in self.contigs.values()]))

	def uncalled_bp(self):
		return(sum([c.seq.count('N') for c in self.contigs.values()]))

class Assembly:
	""" Base class for assembly sequences """
	def __init__(self, args):
		self.name = os.path.basename(args['assembly'])
		self.contigs = self.init_contigs(args)
		self.count_contigs = len(self.contigs.keys())
		self.length = self.length()
		self.uncalled = self.uncalled_bp()
		self.length_trim = self.length - self.uncalled
				
	def init_contigs(self, args):
		contigs = {}
		for rec in Bio.SeqIO.parse(args['assembly'], 'fasta'):
			contigs[rec.id] = Contig(rec)
		return contigs

	def length(self):
		return(sum([contig.length for contig in self.contigs.values()]))

	def uncalled_bp(self):
		return(sum([c.seq.count('N') for c in self.contigs.values()]))

class Contig:
	""" Base class for contig sequences """
	def __init__(self, rec):
		self.name = rec.id
		self.seq = str(rec.seq).upper()
		self.length = len(self.seq)
		self.chunks = []
		self.chunk_ids = set([])

##
## Functions
##

def parse_arguments():
	""" Parse commandline arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="""Evaluate a metagenomic assembly by mapping to a reference genome""")
	parser.add_argument('--assembly', type=str, required=True,
						help="Path to input multi-FASTA file containing one or more contigs")
	parser.add_argument('--genome', type=str, required=False,
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
	parser.add_argument('--pid', type=float, default=99,
						help="Minimum percent identity mapping threshold")
	parser.add_argument('--aln', type=int, default=200,
						help="Minimum alignment length mapping threshold")
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
	cmd += '--fasta %s ' % args['assembly']
	cmd += '--chunk_size %s | ' % args['chunk_size']
	cmd += '%s ' % args['blastn']
	cmd += '-query /dev/stdin '
	cmd += '-subject %s ' % args['genome']
	cmd += '-outfmt 6 '
	cmd += '-max_target_seqs 1 '
	cmd += '-num_threads %s ' % args['threads']
	blastout = run_shell_command(cmd)
	write_m8(args, blastout)

def parse_blast(m8_path):
	""" Yield formatted record from BLAST m8 file """
	formats = {
		0:str, 1:str, 2:float, 3:int, 4:float,
		5:int, 6:int, 7:int, 8:int,
		9:int, 10:float, 11:float
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

def store_alignments(genome, assembly):
	""" Store alignments """
	inpath = os.path.join(args['out'], 'alignments.m8')
	for m8 in parse_blast(inpath):
		assembly_contig = m8['query_id'].rsplit('_', 1)[0]
		chunk_id = m8['query_id'].rsplit('_', 1)[1]
		genome_contig = m8['target_id']
		if m8['query_id'] == 'accn|NC_009513_115':
			print m8
		if m8['aln'] < args['aln']:
			continue
		elif m8['pid'] < args['pid']:
			continue
		elif chunk_id in assembly.contigs[assembly_contig].chunk_ids:
			continue
		else:
			assembly.contigs[assembly_contig].chunk_ids.add(chunk_id)
			assembly.contigs[assembly_contig].chunks.append(m8)
			genome.contigs[genome_contig].chunks.append(m8)

def compute_assembly_stats(assembly):
	""" Compute total alignment length and precision (i.e. chimericity) of entire assembly """
	for contig in assembly.contigs.values():
		contig.aln = sum([c['aln'] for c in contig.chunks])
		contig.ppv = float(contig.aln)/contig.length
	assembly.aln = sum([c.aln for c in assembly.contigs.values()])
	assembly.ppv =  float(assembly.aln)/assembly.length
	assembly.ppv_trim = assembly.aln/float(assembly.length_trim)

def fetch_aln_coords(contig):
	""" Get list of alignment coordinates for each alignment to each contig """
	coords = []
	for c in contig.chunks:
		start = min(c['tstart'], c['tend'])
		end = max(c['tstart'], c['tend'])
		coords.append([start, end])
	return(sorted(coords))

def merge_aln_coords(sorted_coords):
	""" Merge overlapping alignment coordinates """
	if len(sorted_coords) == 0: return []
	merged_coords = [sorted_coords[0]]
	for start, end in sorted_coords[1:]:
		if start <= merged_coords[-1][-1]:
			merged_coords[-1][-1] = end
		else:
			merged_coords.append([start, end])
	return merged_coords

def compute_genome_stats(contig):
	""" Compute number of covered positions on genome """
	for contig in genome.contigs.values():
		sorted_coords = fetch_aln_coords(contig)
		merged_coords = merge_aln_coords(sorted_coords)
		contig.aln = sum([(e - s + 1) for s, e in merged_coords])
		contig.tpr = contig.aln/float(contig.length)
	genome.aln = sum([contig.aln for contig in genome.contigs.values()])
	genome.tpr = genome.aln/float(genome.length)
	genome.tpr_trim = genome.aln/float(genome.length_trim)

def write_report(genome, assembly):
	outfile = open(os.path.join(args['out'], 'summary.txt'), 'w')
	outfile.write('Genome & Assembly Statistics\n')
	outfile.write('%s\t%s\n' % ('genome.name', genome.name))
	outfile.write('%s\t%s\n' % ('genome.count_contigs', genome.count_contigs))
	outfile.write('%s\t%s\n' % ('genome.length', genome.length))
	outfile.write('%s\t%s\n' % ('genome.length_trim', genome.length_trim))
	outfile.write('%s\t%s\n' % ('genome.aln', genome.aln))
	outfile.write('%s\t%s\n' % ('genome.tpr', genome.tpr))
	outfile.write('%s\t%s\n' % ('genome.tpr_trim', genome.tpr_trim))
	outfile.write('\n')
	outfile.write('%s\t%s\n' % ('assembly.name', assembly.name))
	outfile.write('%s\t%s\n' % ('assembly.count_contigs', assembly.count_contigs))
	outfile.write('%s\t%s\n' % ('assembly.length', assembly.length))
	outfile.write('%s\t%s\n' % ('assembly.length_trim', assembly.length_trim))
	outfile.write('%s\t%s\n' % ('assembly.aln', assembly.aln))
	outfile.write('%s\t%s\n' % ('assembly.ppv', assembly.ppv))
	outfile.write('%s\t%s\n' % ('assembly.ppv_trim', assembly.ppv_trim))

##
## Main
##

if __name__ == "__main__":

	args = parse_arguments()
	genome = Genome(args)
	assembly = Assembly(args)

	print("Aligning assembly to reference genome...")
	run_blast(args)
	
	print("Storing results...")
	store_alignments(genome, assembly)
	
	print("Evaluating assembly...")
	compute_assembly_stats(assembly)
	compute_genome_stats(genome)

	print("Writing report...")
	write_report(genome, assembly)



	

