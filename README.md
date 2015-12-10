# Assembly Evaluator  
Evaluate the completedness and precision of a metagenomic assembly by mapping contigs to a reference genome.  
Ideally, reads from reference genome should present in metagenome used to build assembly.

## Dependencies
* Python >=2.6
* BioPython

## Installation
clone software: `git clone https://github.com/snayfach/AssemblyEvaluator`  

## Usage
**usage: evaluate_assembly.py [-h] --assembly FASTA [--genome FASTA] --out DIR  
                            [--threads INT] [--chunk_size INT] [--pid FLOAT]  
							[--aln INT] [--word_size INT]**

Evaluate a metagenomic assembly by mapping to a reference genome

optional arguments:  
 **-h, --help**        show this help message and exit  
  **--assembly** FASTA  Path to input multi-FASTA file containing one or more contigs (default: None)  
  **--genome** FASTA    Path to reference genome FASTA file (default: None)  
  **--out** DIR         Output directory: <out>/alignments.m8, <out>/summary.txt (default: None)  
  **--threads** INT     Threads to use for BLAST (default: 1)  
  **--chunk_size** INT  Chunk size (in bp) for splitting up queries (default:1000)  
  **--pid** FLOAT       Minimum percent identity mapping threshold (default: 99)  
  **--aln** INT         Minimum alignment length mapping threshold (default: 200)  
  **--word_size** INT   Word size for BLAST (default: 50)

## Example
src/evaluate_assembly.py --assembly example/assembly.fa --genome example/genome.fa --out example

## Output file definitions

**genome.name:** file name for reference genome  
**genome.count_contigs:** number of sequences in reference genome  
**genome.length:** total length of reference genome minus ambiguous bases (N)  
**genome.aln:** number of bases in reference genome covered by assembly  
**genome.tpr:** proportion of bases in reference genome covered by assembly

**assembly.name:** file name for metagenomic assembly  
**assembly.count_contigs:** number of contigs in assembly  
**assembly.count_mapped:** number of contigs in assembly that mapped to reference genome  
**assembly.length:** total length of assembly minus ambiguous bases (N)  
**assembly.aln:** total number of bp aligned to reference genome  
**assembly.n50:** N50 of alignment lengths  
**assembly.ppv:** proportion of bp aligned to reference genome for mapped contigs
