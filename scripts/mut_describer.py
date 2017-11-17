#! /usr/local/bin/python2.7 

from collections import Counter, defaultdict
from annotation import mutation, annotation_table
import sys, vcf, pickle, glob
from itertools import permutations
from Bio import SeqIO
from itertools import product, permutations
import gzip

annotation_table_file = sys.argv[1]
reference_fasta = sys.argv[2]
population_vcf_file = sys.argv[3]
mut_quality=0
num_mutated=3
num_unmutated=5
mut_sample_count = defaultdict(int)
h = 0
ref_dict = SeqIO.to_dict(SeqIO.parse(open(reference_fasta, 'r'), 'fasta'))
#fs = glob.glob("CC*/*haplo.indels.txt")
mut_file = sys.argv[4]
vcf_file = sys.argv[5]

print mut_file, mut_file[-2:]


if mut_file[-2:] == "gz":
	f = gzip.open(mut_file,'rb')
else:f = open(mut_file)
#	break_count = 0
line_count = 0
for l in f:
	if line_count == 0 and l[0] == '[': 
		line_count+=1
		continue
	# if break_count > 10: break
	# break_count +=1
	c, p = l.split()[0],int(l.split()[1])
	mut = mutation.mutation(c, p)
	#sys.stderr.write(l)
	mut.ref, mut.alt, mut.qual, mut.MQ, mut.call_method = l.split()[2],  l.split()[3],  float(l.split()[4]), float(l.split()[5]), l.split()[8]
	samples =  mut.samples(vcf_file)
	mutation_event = mut.mutation(vcf_file, mut_quality=mut_quality)
	mut_type = mut.mutation_type(vcf_file)
	sample = mut.mutant_sample(vcf_file, mut_quality=mut_quality)
	mut_sample_count[sample] +=1 #counts mutations per sample over the whole genome
	num_HQ_calls = len(mut.hq_GTs(vcf_file)) #looks for HQ calls - default GQ 20 
	purity = mut.purity(vcf_file)
	if  len(mut.purity(vcf_file)) > 0: min_purity =  min(mut.purity(vcf_file))
	else: min_purity = "n/a"
	het_count = mut.het_count(vcf_file)
	number_mutants = mut.number_mutants(vcf_file)
	number_genotypes =  mut.number_genotypes(vcf_file)
	GTs = mut.GTs(vcf_file)
	ADs = mut.ADs(vcf_file)
	DPs = mut.DPs(vcf_file)
	GQs = mut.GQs(vcf_file)
	gc20 = mut.gc_content(ref_dict, 10)
	gc2000 = mut.gc_content(ref_dict, 1000)
	pi1000 = mut.pi(population_vcf_file, 500, ref_dict)
	annotation_position = mut.annotation_position(annotation_table_file)
	if annotation_position.CDS == '1' and mut_type == 'SNP':
		mut.nonsyn_syn = mut.nonsynonymous(annotation_table_file, vcf_file)
		nonsyn_syn = mut.nonsynonymous(annotation_table_file, vcf_file)
	else: 
		nonsyn_syn = ('n/a','n/a','n/a')
		mut.nonsyn_syn = ('n/a','n/a','n/a')
	nonsense = mut.nonsense(vcf_file, annotation_table_file)
	#header
	annotation_header = "\t".join(['chromosome', \
	'position', \
	'reference_base', \
	'genic', \
	'exonic', \
	'intronic', \
	'intergenic', \
	'utr5', \
	'utr3', \
	'fold0', \
	'fold4', \
	'fold2', \
	'fold3', \
	'CDS ', \
	'mRNA', \
	'rRNA', \
	'tRNA', \
	'feature_names', \
	'feature_types', \
	'feature_ID', \
	'cds_position', \
	'strand', \
	'frame', \
	'codon', \
	'aa', \
	'degen', \
	'FPKM', \
	'rho', \
	'FAIRE', \
	'recombination' \
	])
	#############################################
	header = "\t".join(str(s) for s in \
	[
	'chromosome', \
	'position', \
	'ref', \
	'alt', \
	'qual', \
	'MQ', \
	'mutation', \
	'type', \
	'mutant_sample', \
	'call_method', \
	'HQ_calls', \
	'min_purity', \
	'het_count', \
	'number_mutants', \
	'number_genotypes', \
	'gc20', \
	'gc2000', \
	'pi1000', \
	annotation_header, \
	'Nonsyn_V_Syn\tanc_codon\tmut_codon', \
	'nonsense', \
	"\t".join([str(j)+"_GT" for j in samples]), \
	"\t".join([str(j)+"_AD" for j in samples]), \
	"\t".join([str(j)+"_DP" for j in samples]), \
	"\t".join([str(j)+"_GQ" for j in samples]) \
	])
	if h == 0:
		print header
		h +=1
	#output_data
	sys.stdout.write("\t".join(str(s) for s in \
		[
		mut.chromosome, \
		mut.position, \
		mut.ref, \
		mut.alt, \
		mut.qual, \
		mut.MQ, \
		">".join(mutation_event), \
		mut_type, \
		sample, \
		mut.call_method, \
		num_HQ_calls, \
		min_purity, \
		het_count, \
		number_mutants, \
		number_genotypes, \
		gc20, \
		gc2000, \
		pi1000, \
		annotation_position.output_line(), \
		"\t".join([str(j) for j in nonsyn_syn]), \
		nonsense, \
		"\t".join([str(j) for j in GTs]), \
		"\t".join([str(j) for j in ADs]), \
		"\t".join([str(j) for j in DPs]), \
		"\t".join([str(j) for j in GQs]) \
		]) + "\n")


