#!/usr/bin/python

import os, argparse, gzip, sys, subprocess, Bio.SeqIO, shutil

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
		
	def init_contigs(self, args):
		contigs = {}
		for rec in Bio.SeqIO.parse(args['genome'], 'fasta'):
			contigs[rec.id] = Contig(rec)
		return contigs

	def length(self):
		total_bp = (sum([c.length for c in self.contigs.values()]))
		ambig_bp = uncalled_bp(self.contigs)
		return(total_bp - ambig_bp)

class Assembly:
	""" Base class for assembly sequences """
	def __init__(self, args):
		self.name = os.path.basename(args['assembly'])
		self.contigs = self.init_contigs(args)
		self.count_contigs = len(self.contigs.keys())
		self.length = self.length()
				
	def init_contigs(self, args):
		contigs = {}
		for rec in Bio.SeqIO.parse(args['assembly'], 'fasta'):
			contigs[rec.id] = Contig(rec)
		return contigs

	def length(self):
		total_bp = (sum([c.length for c in self.contigs.values()]))
		ambig_bp = uncalled_bp(self.contigs)
		return(total_bp - ambig_bp)

class Contig:
	""" Base class for contig sequences """
	def __init__(self, rec):
		self.name = rec.id
		self.seq = str(rec.seq).upper()
		self.length = len(self.seq)
		self.m8 = []
		self.chunk_ids = set([])
		self.mapped = False

##
## Functions
##

def uncalled_bp(contigs):
	""" Count the number of ambiguous base calls across contigs """
	return(sum([c.seq.count('N') for c in contigs.values()]))

def parse_arguments():
	""" Parse commandline arguments """
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="""Evaluate a metagenomic assembly by mapping to a reference genome""")
	parser.add_argument('--assembly', type=str, required=True, metavar='FASTA',
						help="Path to input multi-FASTA file containing one or more contigs")
	parser.add_argument('--genome', type=str, required=False, metavar='FASTA',
						help="Path to reference genome FASTA file")
	parser.add_argument('--out', dest='out', type=str, required=True, metavar='DIR',
						help="""
						Output directory: <out>/alignments.m8, <out>/summary.txt
						""")
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
						help="Threads to use for BLAST")
	parser.add_argument('--chunk_size', type=int, default=1000, metavar='INT',
						help="Chunk size (in bp) for splitting up queries")
	parser.add_argument('--pid', type=float, default=99, metavar='FLOAT',
						help="Minimum percent identity mapping threshold")
	parser.add_argument('--aln', type=int, default=200, metavar='INT',
						help="Minimum alignment length mapping threshold")
	parser.add_argument('--word_size', type=int, default=50, metavar='INT',
						help="Word size for BLAST")
	args = vars(parser.parse_args())
	if not os.path.isdir(args['out']): os.mkdir(args['out'])
	add_binaries(args)
	return args

def add_binaries(args):
	""" Add paths to binaries """
	script_path = os.path.abspath(__file__)
	script_dir = os.path.dirname(script_path)
	main_dir = os.path.dirname(script_dir)
	args['hs-blastn'] = '%s/bin/hs-blastn' % main_dir
	args['blastn'] = '%s/bin/blastn' % main_dir
	args['makeblastdb'] = '%s/bin/makeblastdb' % main_dir
	args['stream_chunks'] = '%s/stream_chunks.py' % script_dir
	args['db'] = os.path.join(args['out'], os.path.basename(args['genome']))

def run_shell_command(command):
	""" Capture stdout, stderr. Check unix exit code and exit if non-zero """
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s" % (command, err)
		sys.exit(err_message)
	else:
		return out

def build_blastdb(args):
	""" Build HS-BLASTN database """
	if not os.path.isfile(args['db']):
		shutil.copy(args['genome'], args['db'])
	cmd = '%s index %s ' % (args['hs-blastn'], args['db'])
	out = run_shell_command(cmd)

def run_blast(args):
	""" Run HS-BLASTN """
	cmd = '%s ' % args['stream_chunks']
	cmd += '--fasta %s ' % args['assembly']
	cmd += '--chunk_size %s | ' % args['chunk_size']
	cmd += '%s align ' % args['hs-blastn']
	cmd += '-query /dev/stdin '
	cmd += '-db %s ' % args['db']
	cmd += '-outfmt 6 '
	cmd += '-num_alignments 1 '
	cmd += '-word_size %s ' % args['word_size']
	cmd += '-num_threads %s ' % args['threads']
	blastout = run_shell_command(cmd)
	return(blastout)

def write_m8(args, blastout):
	""" Write tab delimited m8 file to disk """
	outfile = open(os.path.join(args['out'], 'alignments.m8'), 'w')
	outfile.write(blastout)
	outfile.close()

def blast_cleanup(args):
	""" Clean blast temporary files """
	basename = os.path.join(args['out'], os.path.basename(args['genome']))
	exts = ['', '.bwt', '.header', '.sa', '.sequence']
	for ext in exts:
		os.remove(basename+ext)

def align_seqs(args):
	""" Run HS-BLASTN """
	build_blastdb(args)
	blastout = run_blast(args)
	write_m8(args, blastout)
	blast_cleanup(args)

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

def store_alignments(genome, assembly):
	""" Store alignments """
	inpath = os.path.join(args['out'], 'alignments.m8')
	for m8 in parse_blast(inpath):
		assembly_contig = m8['query_id'].rsplit('_', 1)[0]
		chunk_id = m8['query_id'].rsplit('_', 1)[1]
		genome_contig = m8['target_id']
		if m8['aln'] < args['aln']:
			continue
		elif m8['pid'] < args['pid']:
			continue
		elif chunk_id in assembly.contigs[assembly_contig].chunk_ids:
			continue
		else:
			assembly.contigs[assembly_contig].mapped = True
			assembly.contigs[assembly_contig].chunk_ids.add(chunk_id)
			assembly.contigs[assembly_contig].m8.append(m8)
			genome.contigs[genome_contig].m8.append(m8)

def assembly_ppv(assembly):
	""" Use mapped contigs to compute assembly precision """
	length = 0
	aln = 0
	for contig in assembly.contigs.values():
		if contig.mapped:
			length += contig.length
			aln += contig.aln
	return float(aln)/length

def compute_assembly_stats(assembly):
	""" Compute total alignment length and precision (i.e. chimericity) of entire assembly """
	for contig in assembly.contigs.values():
		contig.aln = sum([c['aln'] for c in contig.m8])
		contig.ppv = float(contig.aln)/contig.length
	assembly.mapped_length = sum([c.length if c.mapped else 0 for c in assembly.contigs.values()])
	assembly.aln = sum([c.aln for c in assembly.contigs.values()])
	assembly.ppv =  assembly_ppv(assembly)
	assembly.n50 = assembly_n50(assembly)
	assembly.count_mapped = sum([c.mapped for c in assembly.contigs.values()])

def fetch_aln_coords(contig):
	""" Get list of alignment coordinates for each alignment to each contig """
	coords = []
	for c in contig.m8:
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
		contig.alns = [(e - s + 1) for s, e in merged_coords]
		contig.aln = sum(contig.alns)
		contig.tpr = contig.aln/float(contig.length)
	genome.aln = sum([contig.aln for contig in genome.contigs.values()])
	genome.tpr = genome.aln/float(genome.length)

def assembly_n50(assembly):
	""" Compute n50 based on contig alignment lengths """
	x = 0
	sorted_alns = sorted([c.aln for c in assembly.contigs.values()], reverse=True)
	for a in sorted_alns:
		if x >= assembly.length/2.0: return(a)
		else: x += a
	return a

def write_asm_report(genome, assembly):
	""" Report alignment results for genome and assembly """
	outfile = open(os.path.join(args['out'], 'asm.summary.txt'), 'w')
	outfile.write('%s\t%s\n' % ('genome.name', genome.name))
	outfile.write('%s\t%s\n' % ('genome.count_contigs', genome.count_contigs))
	outfile.write('%s\t%s\n' % ('genome.length', genome.length))
	outfile.write('%s\t%s\n' % ('genome.aln', genome.aln))
	outfile.write('%s\t%s\n' % ('genome.tpr', genome.tpr))
	outfile.write('\n')
	outfile.write('%s\t%s\n' % ('assembly.name', assembly.name))
	outfile.write('%s\t%s\n' % ('assembly.count_contigs', assembly.count_contigs))
	outfile.write('%s\t%s\n' % ('assembly.count_mapped', assembly.count_mapped))
	outfile.write('%s\t%s\n' % ('assembly.total_length', assembly.length))
	outfile.write('%s\t%s\n' % ('assembly.mapped_length', assembly.mapped_length))
	outfile.write('%s\t%s\n' % ('assembly.aln', assembly.aln))
	outfile.write('%s\t%s\n' % ('assembly.n50', assembly.n50))
	outfile.write('%s\t%s\n' % ('assembly.ppv', assembly.ppv))

def write_contig_report(assembly):
	""" Report alignment results for assembled contigs """
	outfile = open(os.path.join(args['out'], 'contig.summary.txt'), 'w')
	fields = ['contig.name', 'contig.length', 'contig.mapped', 'contig.aln', 'contig.ppv']
	outfile.write('\t'.join([str(_) for _ in fields])+'\n')
	for contig in assembly.contigs.values():
		values = [contig.name, contig.length, contig.mapped, contig.aln, contig.ppv]
		outfile.write('\t'.join([str(_) for _ in values])+'\n')

##
## Main
##

if __name__ == "__main__":

	args = parse_arguments()
	genome = Genome(args)
	assembly = Assembly(args)

	print("Aligning assembly to reference genome...")
	align_seqs(args)
	
	print("Storing results...")
	store_alignments(genome, assembly)
	
	print("Evaluating assembly...")
	compute_assembly_stats(assembly)
	compute_genome_stats(genome)

	print("Writing reports...")
	write_asm_report(genome, assembly)
	write_contig_report(assembly)


	

