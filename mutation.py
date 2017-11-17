#! /usr/bin/python
from string import translate, maketrans
import sys, vcf, re
from Bio import SeqIO
from operator import itemgetter, attrgetter
from pysam import Tabixfile
from ness_vcf import ness_vcf #pi_from_AF, consensus_sequence, indel, AF2pi, AF1
from itertools import permutations
from collections import Counter, OrderedDict
from annotation import annotation_table
############################################################################################################################################
def reverse_complement(sequence):
	return str(sequence)[::-1].translate(maketrans('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))
############################################################################################################################################
class mut_line(object):
	def __init__(self, line):
		self.line = line.strip()
		self.chromosome, position, self.ref, self.alt, qual, GTs, GQs, ADs = line.split("\t")
		self.position = int(position)
		self.qual = float(qual)
		exec "GTs=[%s]" %( ",".join(["'"+i+"'" for i in GTs.split(",")]))
		exec "GQs=[%s]" %(GQs)
		exec "ADs=[%s]" %(ADs)
		self.GTs = GTs
		self.GQs = GQs
		self.ADs = ADs
		self.purities = False
	def purity(self):
		if self.purities: return self.purities
		self.purities = [float(max(i))/sum(i) for i in self.ADs if i!= None and sum(i) > 0]
		return self.purities
	def het_count(self):
		self.het_count = 0
		for i in self.GTs:
			if i != 'None' and i.split("/")[0] != i.split("/")[1]:
				self.het_count +=1
		return self.het_count
	def fishy(self):
		result =  [0,0]
		g = [i for i in self.GTs if (i and i != 'None')]
		if len(Counter(g)) > 2:
			result[0] = 1
			#return True
		if Counter(g).most_common(2)[1][1] > 1: 
			#print "Multiple individuals with mutant genotype", g
			result[1] = 1
			#return True
		return result
		# else:
		# 	return "fine"
		# 	#return False


############################################################################################################################################
############################################################################################################################################

class mutation(object):
	""" 
	To fully describe a mutation it must have a number of files
	it must have its accompanying VCF
	it could have its accompanying reference_fasta
	it must have its accompanying annotation_table
	its initialization could only really be chromosome and position
	Things to collect
	genome_position_object
	gene density fraction of sites in window that are "genic/exonic/CDS"
	distance to exon go up and down until you hit an exon
	distance to origin of replication # no idea
	genic annotation
		coding, intronic, UTR5, UTR3, non-synonymous, synonymous, exon number, position in gene
	expression
	diversity short window, large window
	Neutral diversity?
	Divergence where possible
	clustering
	transition/transversion
	Genotype
	sequencing stats
	coverage
	"""
	def __init__(self, chromosome, position):
		self.chromosome = chromosome
		self.position = position
		#things that could be defined 
		self.strain = False #this can just be set, ie CC1373, 2937 etc
		self.annotation_table_file = False #this is a file name for the genome annotation table
		self.vcf_file = False
		self._mutant_genotype = False 
		self.GQ = False # I don't really know what I was thinking here
		self._mutation_type = False #indel or snp
		self._mutation = False #(vcf_record.REF, vcf_record.ALT) #note this should account for the ancestral base and derived base
		self._mutant_genotype = False
		self._mutant_sample=False
		self._annotation_position = False #booyah
		self._vcf_record = False
		self._vcf_line = False
		self._gc_content = False
		self._pi= False
		self._gene_density = False
		self._hq_GTs = False
		self._hq_GQs = False
		self._hq_samples = False
		self._hq_ADs = False
		self._hq_DPs = False
		self._ADs = False
		self._DPs = False
		self._GTs = False
		self._GQs = False
		self._DPs = False
		self._samples = False
		self._purity = False
		self._het_count = False
		self._number_mutants = False
		self._number_genotypes = False
		self._is_variant = False
		self._upstream_region = False
	#
	def annotation_position(self, annotation_table_file):
		"""Note for this function to work the annotation table has to be tabix'ed"""
		if self._annotation_position: return self._annotation_position
		else:
			self._annotation_position = annotation_table.annotation_table_position(self.chromosome, self.position, annotation_table_file)
			return self._annotation_position

	def vcf_line(self, vcf_file):
		#the problem with this is that there may be two VCF lines with this chromosome position
		if self._vcf_line: return self._vcf_line
		else:
			v = Tabixfile(vcf_file)
			self._vcf_line = v.fetch(reference=self.chromosome, start=self.position-1, end=self.position).next()
			v.close()
			return self._vcf_line
	
	def vcf_record(self, vcf_file, num_mutated=1, num_unmutated=1, mut_quality=0):
		if self._vcf_record: return self._vcf_record
		else:
			defined = False
			vr = vcf.Reader(filename=vcf_file)
			for v in vr.fetch(self.chromosome, self.position-1, self.position):
				if v.is_snp or v.is_indel and v.POS == self.position:
					hq_GTs = [s['GT'] for s in v.samples if s['GQ']>mut_quality ]  #self.hq_GTs(vcf_file)
					hq_GQs = [s['GQ'] for s in v.samples if s['GQ']>mut_quality ] #self.hq_GQs(vcf_file)
					for p in [q for q in permutations(range(len(v.ALT)+1))]: #this basically makes all pairs of all the possible alleles at a site 0,1 0,2 1,2 2,1
						mut = "%i/%i" %(p[0], p[0])
						wt = "%i/%i" %(p[1], p[1])
						if num_mutated>=hq_GTs.count(mut)>0 and hq_GTs.count(wt) >= num_unmutated:
							self._vcf_record = v
							defined = True
							return self._vcf_record
			if not defined:
				#still can't find it
				for v in vr.fetch(self.chromosome, self.position, self.position):
					if v.is_snp or indel(v) and v.POS == self.position:
						self._vcf_record = v
						return self._vcf_record
					else: print "BOOM"
	
	def mutant_genotype(self, vcf_file, num_mutated=1, num_unmutated=1, mut_quality=20):
		if self._mutant_genotype: return self._mutant_genotype
		else:
			vr = self.vcf_record(vcf_file, mut_quality=mut_quality) 
			hq_GTs = self.hq_GTs(vcf_file, mut_quality=mut_quality)
			hq_GQs = self.hq_GQs(vcf_file, mut_quality=mut_quality)
			for p in [q for q in permutations(range(len(vr.ALT)+1))]: #this basically makes all pairs of all the possible alleles at a site 0,1 0,2 1,2 2,1
				mut = "%i/%i" %(p[0], p[0])
				wt = "%i/%i" %(p[1], p[1])
				if num_mutated>=hq_GTs.count(mut)>0 and hq_GTs.count(wt) >= num_unmutated:
					if self.GTs(vcf_file).count(mut)>self.GTs(vcf_file).count(wt):
						#something has gone wrong. Often there are only 2 HQ calls and the wrong one has been chosen as the mutant
						#I can normally rectify this by just reversing the mut and wt genotypes
						#however, in some cases its just not a mutation.
						#if after reversing the count of the mutation is still over 1 its probably not a mutation.
						self._mutant_genotype = wt
					else:
						self._mutant_genotype = mut
					return self._mutant_genotype
	
	def hq_GTs(self, vcf_file, mut_quality=20):
		if self._hq_GTs:return self._hq_GTs
		else:
			self._hq_GTs = [s['GT'] for s in self.vcf_record(vcf_file).samples if s['GQ']>mut_quality ]
			return self._hq_GTs
	
	def hq_GQs(self, vcf_file, mut_quality=20):
		if self._hq_GQs:return self._hq_GQs
		else:
			self._hq_GQs = [s['GQ'] for s in self.vcf_record(vcf_file).samples if s['GQ']>mut_quality ]
			return self._hq_GQs
	
	def hq_samples(self, vcf_file, mut_quality=20):
		if self._hq_samples:return self._hq_samples
		else:
			self._hq_samples = [s.sample for s in self.vcf_record(vcf_file).samples if s['GQ']>mut_quality ]
			return self._hq_samples
	
	def hq_ADs(self, vcf_file, mut_quality=20):
		if self._hq_ADs:return self._hq_ADs
		else:
			self._hq_ADs = [s['AD'] for s in self.vcf_record(vcf_file).samples if s['GQ']>mut_quality ]
			return self._hq_ADs
	
	def hq_DPs(self, vcf_file, mut_quality=20):
		if self._hq_DPs:return self._hq_DPs
		else:
			self._hq_DPs = [s['DP'] for s in self.vcf_record(vcf_file).samples if s['GQ']>mut_quality ]
			return self._hq_DPs
	
	def GTs(self, vcf_file):
		if self._GTs:return self._GTs
		else:
			self._GTs = [s['GT'] for s in self.vcf_record(vcf_file).samples]
			return self._GTs
	
	def GQs(self, vcf_file):
		if self._GQs:return self._GQs
		else:
			self._GQs = [s['GQ'] for s in self.vcf_record(vcf_file).samples]
			return self._GQs
	
	def samples(self, vcf_file):
		if self._samples:return self._samples
		else:
			self._samples = [s.sample for s in self.vcf_record(vcf_file).samples]
			return self._samples
	
	def ADs(self, vcf_file):
		if self._ADs:return self._ADs
		else:
			self._ADs = [s['AD'] for s in self.vcf_record(vcf_file).samples]
			return self._ADs
	
	def DPs(self, vcf_file):
		if self._DPs:return self._DPs
		else:
			self._DPs = [s['DP'] for s in self.vcf_record(vcf_file).samples]
			return self._DPs
	
	def purity(self, vcf_file):
		if self._purity: return self._purity
		else:
			self._purity = [float(max(i))/sum(i) for i in self.ADs(vcf_file) if i!= None and sum(i) > 0]
			return self._purity
	
	def het_count(self, vcf_file):
		if self._het_count: return self._het_count
		else:
			self._het_count = 0
			for i in self.GTs(vcf_file):
				if i != None and i.split("/")[0] != i.split("/")[1]:
					self._het_count +=1
			return self._het_count
	
	def mutation(self, vcf_file, mut_quality=20, num_mutated=1, num_unmutated=1):
		# should I make this take into account the consensus base????) -- yes
		if self._mutation: return self._mutation
		vcf_record = self.vcf_record(vcf_file)
		hq_GTs = self.hq_GTs(vcf_file, mut_quality= mut_quality)
		hq_GQs = self.hq_GQs(vcf_file, mut_quality= mut_quality) 
		for p in [i for i in permutations(range(len(vcf_record.ALT)+1))]: #this basically makes all pairs of all the possible alleles at a site 0,1 0,2 1,2 2,1
			mut = "%i/%i" %(p[0], p[0])
			wt = "%i/%i" %(p[1], p[1])
			if num_mutated>=hq_GTs.count(mut)>0 and hq_GTs.count(wt) >= num_unmutated:
				if p[0] == 0:
					mut_base = vcf_record.REF
				else:
					mut_base=vcf_record.ALT[p[0]-1]
				if p[1] == 0: 
					ancestral_base = vcf_record.REF
				else:
					ancestral_base=vcf_record.ALT[p[1]-1]
				break #this stops once its found a mutation
		if self.GTs(vcf_file).count(mut)>self.GTs(vcf_file).count(wt):
			self._mutation = (str(mut_base), str(ancestral_base))
		else:
			self._mutation = (str(ancestral_base), str(mut_base))	
		return self._mutation
	
	def mutant_sample(self, vcf_file, num_mutated=1, num_unmutated=1, mut_quality=20):
		if self._mutant_sample: return self._mutant_sample
		else:
			hq_GTs = self.hq_GTs(vcf_file, mut_quality= mut_quality) 
			hq_GQs = self.hq_GQs(vcf_file, mut_quality= mut_quality) 
			hq_samples = self.hq_samples(vcf_file, mut_quality= mut_quality)
			for p in [q for q in permutations(range(len(self.vcf_record(vcf_file).ALT)+1))]: #this basically makes all pairs of all the possible alleles at a site 0,1 0,2 1,2 2,1
				mut = "%i/%i" %(p[0], p[0])
				wt = "%i/%i" %(p[1], p[1])
				if num_mutated>=hq_GTs.count(mut)>0 and hq_GTs.count(wt) >= num_unmutated:
					e = 0
					for s in range(len(hq_GTs)):
						if hq_GTs[s] == mut:
							self._mutant_sample = hq_samples[s]
							return self._mutant_sample
	
	def number_mutants(self, vcf_file):
		if self._number_mutants: return self._number_mutants
		else:
			self._number_mutants = self.GTs(vcf_file).count(self.mutant_genotype(vcf_file))
			return self._number_mutants
	
	def number_genotypes(self, vcf_file):
		if self._number_genotypes: return self._number_genotypes
		else:
			gts = [i for i in self.GTs(vcf_file) if i != None]
			self._number_genotypes = len(set(gts)) 
			return self._number_genotypes
	
	def gc_content(self, reference_dict, window_size):
		window_start = max(0,self.position-window_size)
		window_end = min(len(reference_dict[self.chromosome].seq),self.position+window_size+1)
		window_sequence = str(reference_dict[self.chromosome].seq[window_start:window_end])
		gc_content = len(re.findall('[GCgc]', window_sequence))/float(len(re.findall('[ACGTacgt]', window_sequence)))
		return gc_content
	
	def pi(self, vcf_file, window_size, reference_dict, min_called=2, ploidy=1 ):
		window_start = max(0,self.position-window_size)
		window_end = min(len(reference_dict[self.chromosome].seq),self.position+window_size+1)
		vcf_reader = vcf.Reader(filename=vcf_file)
		region = "%s:%i-%i" %(self.chromosome, window_start,window_end)
		pi = ness_vcf.pi_from_AF(vcf_reader, region, min_called=min_called, ploidy=ploidy)[-1]
		return pi
	
	def gene_density(self, window_size, annotation_table_file, reference_dict, feature_type='genic'):
		window_start = max(0,self.position-window_size)
		window_end = min(len(reference_dict[self.chromosome].seq),self.position+window_size+1)
		genic_sites=0
		for genome_position in range(window_start,window_end): #this could potentially run off the end of a chromosome
			#this looks wrong.
			p =  annotation_table.annotation_table_position(self.chromosome, genome_position, annotation_table_file )
			genic_sites +=  int(p.genic)
		# more concise but no faster
		# genic_sites = sum([int(annotation_table_position(self.chromosome, genome_position,annotation_table_file).genic) \
		# 					for genome_position in range(max(1,(self.position-window_size)),self.position+window_size)])
		self._gene_density = genic_sites/(window_size*2.0)
		return self._gene_density
	
	def mutation_type(self, vcf_file):
		if self._mutation_type: return self._mutation_type
		mutation_event = self.mutation(vcf_file)
		if len(mutation_event[0]) > 1 or len(mutation_event[1]) > 1:
			self._mutation_type = 'indel'
		else: self._mutation_type = 'SNP'
		return self._mutation_type 
		# if self.vcf_record(vcf_file).is_snp:
		# 	self._mutation_type = 'SNP'
		# 	return self._mutation_type
		# elif indel(self.vcf_record(vcf_file)):
		# 	self._mutation_type = 'indel'
		# 	return self._mutation_type
	
	def nonsynonymous(self, annotation_table_file, vcf_file):
		"""
		determine if a particular mutation changes the encoded protein
		must take into account the ancestral sequence and the codon
			to get the ancestral sequence you need to retrieve the positions of the codon
			once you have the positions get the consensus bases for that position to assemble the codon
			then take the mutation and determine whether it causes a non-syn mutation
		"""
		gen_code_dict = { 'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', \
			'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S', \
			'tat': 'Y', 'tac': 'Y', 'taa': 'X', 'tag': 'X', \
			'tgt': 'C', 'tgc': 'C', 'tga': 'X', 'tgg': 'W', \
			'ctt': 'L', 'ctc': 'L', 'cta': 'L', 'ctg': 'L', \
			'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', \
			'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q', \
			'cgt': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R', \
			'att': 'I', 'atc': 'I', 'ata': 'I', 'atg': 'M', \
			'act': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T', \
			'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K', \
			'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R', \
			'gtt': 'V', 'gtc': 'V', 'gta': 'V', 'gtg': 'V', \
			'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', \
			'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E', \
			'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'}
		def codon_positions(annotation_position, annotation_table_file):
			a=annotation_position
			def get_it(chromosome, start_position, direction, frame, annotation_table_file):
				"""
				This helper function returns the genome position 
				which matches 'frame' by moving up (-) or down (+) the genome 
				in the annotation table in the direction indicated
				"""
				p = annotation_table.annotation_table_position(chromosome, int(start_position), annotation_table_file)
				while p.frame[0] != frame:
					if direction == "-": new_position = int(p.position) - 1
					elif direction == "+": new_position = int(p.position) + 1
					p = annotation_table.annotation_table_position(chromosome, new_position, annotation_table_file)
				#now we have the focal position stored as p
				return p
			if a.strand[0] == '-':
				if a.frame[0] == '0':
					p0 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p1 = get_it(a.chromosome, a.position, "-", '1', annotation_table_file)
					p2 = get_it(a.chromosome, a.position, "-", '2', annotation_table_file)
				elif a.frame[0] == '1':
					p1 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p2 = get_it(a.chromosome, a.position, "-", '2', annotation_table_file)
					p0 = get_it(a.chromosome, a.position, "+", '0', annotation_table_file)
				elif a.frame[0] == '2':
					p2 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p1 = get_it(a.chromosome, a.position, "+", '1', annotation_table_file)
					p0 = get_it(a.chromosome, a.position, "+", '0', annotation_table_file)
				return (p2, p1,p0) #ie this is not reversed or complemented
			elif a.strand[0] == "+":
				if a.frame[0] == '0':
					p0 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p1 = get_it(a.chromosome, a.position, "+", '1', annotation_table_file)
					p2 = get_it(a.chromosome, a.position, "+", '2', annotation_table_file)
				elif a.frame[0] == '1':
					p1 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p2 = get_it(a.chromosome, a.position, "+", '2', annotation_table_file)
					p0 = get_it(a.chromosome, a.position, "-", '0', annotation_table_file)
				elif a.frame[0] == '2':
					p2 = annotation_table.annotation_table_position(a.chromosome, a.position, annotation_table_file)
					p1 = get_it(a.chromosome, a.position, "-", '1', annotation_table_file)
					p0 = get_it(a.chromosome, a.position, "-", '0', annotation_table_file)
				return (p0,p1,p2)
			else:
				print "what is this strand:", a.strand
				sys.exit()
		#this is the annotation line of the mutation
		a = self.annotation_position(annotation_table_file)
		positions = codon_positions(a, annotation_table_file)
		chromosome = a.chromosome
		mutation_event = self.mutation(vcf_file)
		if a.strand[0] == "-":
			codon = reverse_complement(''.join([ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[0].position), \
			ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[1].position), \
			ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[2].position)]))
			#swap out the affected base
			mutated_codon = list(codon) #change it to a list so we can 
			mutated_codon[int(a.frame[0])] = reverse_complement(mutation_event[1])
			mutated_codon = ''.join(mutated_codon)
		elif a.strand[0] == "+":
			codon = ''.join([ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[0].position), \
			ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[1].position), \
			ness_vcf.consensus_sequence(vcf_file, self.chromosome, positions[2].position)])
			mutated_codon = list(codon) #change it to a list so we can 
			mutated_codon[int(a.frame[0])] = mutation_event[1]
			mutated_codon = ''.join(mutated_codon)
		if codon.lower() in gen_code_dict.keys():
			if gen_code_dict[codon.lower()] != gen_code_dict[mutated_codon.lower()]:
				return ('nonsynonymous', codon, mutated_codon)
				return True
			elif gen_code_dict[codon.lower()] == gen_code_dict[mutated_codon.lower()]:
				return ('synonymous', codon, mutated_codon)
				return False
		else: 
			return ('NA', codon, mutated_codon)
			return False
	
	def nonsense(self, vcf_file, annotation_table_file):
		"""check for indels in CDS that aren't divisible by three"""
		mutation_event = self.mutation(vcf_file)
		annotation_position = self.annotation_position(annotation_table_file)
		if annotation_position.CDS == "1" and (len(mutation_event[0]) > 1 or len(mutation_event[1]) > 1):
			indel_len = abs(len(mutation_event[0]) - len(mutation_event[1]))
			if indel_len%3 != 0:
				return 'nonsense'
			else: 
				return 'sense'
		else:
			return 'n/a'
	
	def TsTv(self, vcf_file):
		base_dict = {'A':'purine', 'G':'purine', 'C':'pyrimidine', 'T':'pyrimidine'}
		ancestral_base, mut_base = self.mutation(vcf_file)
		if ancestral_base in base_dict and mut_base in base_dict:
			if base_dict[ancestral_base] == base_dict[mut_base]: return 'transition'
			else: return 'transversion'
	def dist2centromere(self):
 		centromere_dict = {'chromosome_1': 835000,
							'chromosome_2': 2470000,
							'chromosome_3': 6570000,
							'chromosome_4': 713634,
							'chromosome_5': 2061000,
							'chromosome_6': 4120000,
							'chromosome_7': 3488000,
							'chromosome_8': 2520000,
							'chromosome_9': 4559000,
							'chromosome_10': 2927000,
							'chromosome_11': 2546000,
							'chromosome_12': 7150000,
							'chromosome_13': 4085000,
							'chromosome_14': 2693000,
							'chromosome_15': 1097941,
							'chromosome_16': 3776000,
							'chromosome_17': 6111000}
 		if self.chromosome in centromere_dict: return abs(centromere_dict[self.chromosome] - self.position)
		else: return None
	
	def is_variant(self, vcf_file):
		"""check if site is variant in population vcf"""
		if self._is_variant: return self._is_variant
		r = vcf.Reader(filename=vcf_file).fetch(self.chromosome, self.position, self.position).next()
		if r.is_snp or indel(r):
			self._is_variant = True
			return self._is_variant
	
	def upstream_region(self, region_size, annotation_table_file):
		ch_data ={'cpDNA': 203828, 'scaffold_30': 52376, 'scaffold_51': 11225, 'scaffold_33': 39192, 'scaffold_53': 2479, 'scaffold_52': 6241, 'scaffold_54': 2277, 'scaffold_32': 42264, 'scaffold_31': 48183, 'chromosome_15': 1922860, 'chromosome_14': 4157777, 'chromosome_17': 7188315, 'chromosome_16': 7783580, 'chromosome_11': 3826814, 'chromosome_10': 6576019, 'chromosome_13': 5206065, 'chromosome_12': 9730733, 'scaffold_37': 24537, 'scaffold_39': 22408, 'scaffold_38': 24437, 'scaffold_19': 219038, 'scaffold_18': 271631, 'scaffold_36': 25399, 'scaffold_35': 32450, 'scaffold_50': 12727, 'scaffold_34': 33576, 'mtMinus': 345555, 'mtDNA': 15758, 'chromosome_5': 3500558, 'chromosome_4': 4091191, 'chromosome_7': 6421821, 'chromosome_6': 9023763, 'chromosome_1': 8033585, 'chromosome_3': 9219486, 'chromosome_2': 9223677, 'chromosome_9': 7956127, 'chromosome_8': 5033832, 'scaffold_46': 16627, 'scaffold_47': 14746, 'scaffold_44': 17736, 'scaffold_45': 16939, 'scaffold_42': 21000, 'scaffold_43': 20974, 'scaffold_40': 22082, 'scaffold_41': 21325, 'scaffold_48': 14165, 'scaffold_49': 13462, 'scaffold_20': 200793, 'scaffold_21': 189560, 'scaffold_22': 163774, 'scaffold_23': 127913, 'scaffold_24': 127161, 'scaffold_25': 102191, 'scaffold_26': 80213, 'scaffold_27': 55320, 'scaffold_28': 55278, 'scaffold_29': 52813}
		if self._upstream_region: return self._upstream_region
			#so basically walk up and down the annotation table and see whether a site in that region is UTR5 or beginning of exon
		if self.annotation_position(annotation_table_file).genic == '0':
			for site in range(max(1,self.position-region_size), self.position):
				a = annotation_table.annotation_table_position(self.chromosome, site, annotation_table_file)
				if a.chromosome == self.chromosome and a.genic == '1' and  (a.strand == "-" or a.utr5 =='1'):
					#we found the start of a a gene
					self._upstream_region = '1'
					return self._upstream_region
			for site in range(self.position, min(self.position+region_size,ch_data[self.chromosome])):
				a = annotation_table.annotation_table_position(self.chromosome, site, annotation_table_file)
				if a.chromosome == self.chromosome and a.genic == '1' and (a.strand == "+" or a.utr5 =='1'):
					self._upstream_region = '1'
					return self._upstream_region
			if self._upstream_region == False: 
				self._upstream_region = '0'
				return self._upstream_region 
		else:
			self._upstream_region = '0'
			return self._upstream_region 

def pi_from_AF_by_site(vcf_reader, region, annotation_table_file, site_category, min_called=2, ploidy=1):
	"""  calculates theta pi for a region of vcf """
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	num_alleles = ploidy*len(vcf_reader.samples)
	def check_site(annotation_table_file, chromosome, position, site_category):
		if eval('annotation_table_position(chromosome, position, annotation_table_file).' + site_category) == '1':return True
		else: return False
	pis =  [AF2pi(AF1(i, min_called)[0],num_alleles) for i in vcf_reader.fetch(chromosome, start_coord, end_coord) if AF1(i, min_called)!=None and check_site(annotation_table_file, chromosome, i.POS, site_category)]
	return [region, len([i for i in pis if i>0]), len(pis), sum(pis)/len(pis)]

