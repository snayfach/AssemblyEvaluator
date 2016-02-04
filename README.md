# Assembly Evaluator  
Evaluate the completedness and precision of a (meta)genomic assembly by mapping contigs to a complete reference genome.  

Use this tool to evaluate an assembly when you know *apriori* that a reference genome is present. For example, if you "spiked-in" simulated reads, or if you are assembling an in-silico simulated metagenome.

This tool will map the contigs back to the reference genome and answer three questions:
* (Completedness) What fraction of the reference is covered by mapped contigs?
* (Precision/Chimericity) Of those mapped contigs, what fraction is aligned to the reference?
* (Quality) How long are the mapped contigs? What is their N50?

Use cases:
* Spike-in short-reads from a reference genome into a complex metagenome and assemble data into contigs. Use AssemblyEvaluator to evaluate the resulting contigs for precision and completedness
* Other cases when assembling data where is a completed reference genome to map to
* You can also use the tool to evaluate a set of genes identified from assembled contigs

## Dependencies
* Python >=2.6
* BioPython

## Installation
clone software: `git clone https://github.com/snayfach/AssemblyEvaluator`  

## Usage
**evaluate_assembly.py [-h] --assembly FASTA --genome FASTA --out DIR  
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

* evaluate assembled contigs versus a reference genome:  
`src/evaluate_assembly.py --assembly example/assembled_contigs.fa --genome example/reference_genome.fa --out example`

* evaluate assembled genes versus a genes from a reference genome:  
`src/evaluate_assembly.py --assembly example/assembled_genes.fa --genome example/reference_genes.fa --out example`


## Output files

### asm.summary.txt
*genome.name:* file name for reference genome  
*genome.count_contigs:* number of sequences in reference genome  
*genome.length:* total length of reference genome minus ambiguous bases (N)  
*genome.aln:* number of bases in reference genome covered by assembly  
*genome.tpr:* proportion of bases in reference genome covered by assembly

*assembly.name:* file name for metagenomic assembly  
*assembly.count_contigs:* number of contigs in assembly  
*assembly.count_mapped:* number of contigs in assembly that mapped to reference genome  
*assembly.total_length:* length of assembly minus ambiguous bases (N)  
*assembly.mapped_length:* length of mapped contigs minus ambiguous bases (N)    
*assembly.aln:* total number of bp aligned to reference genome  
*assembly.n50:* N50 of alignment lengths  
*assembly.ppv:* proportion of bp aligned to reference genome for mapped contigs

The most important variables are **genome.tpr**, which indicates the fraction of reference genome/genes covered by the assembly, and **assembly.ppv**, which indicates the fraction of mapped contig length that is aligned to the genome/genes.

A low **genome.tpr** indicates that most of the reference genome/genes was not recovered by the assembly

A low **assembly.ppv** indicates that a significant proportion of assembled contigs that mapped to the reference were chimeric

### contig.summary.txt
Per-contig report

### alignments.m8
tabular outpur from HS-BLASTN
