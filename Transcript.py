import sys, vcf, re
from Bio import SeqIO
from operator import itemgetter, attrgetter
from itertools import permutations
from collections import Counter, OrderedDict
from annotation import GFF_line


#########################################################################################################################################
def reverse_complement(sequence):
	#python 2 return str(sequence)[::-1].translate(maketrans('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))
	tr = dict(zip('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))
	return "".join([tr[i.upper()] for i in sequence])[::-1]
#########################################################################################################################################

class Output_Line(): 
	#this is a class that stores all the field that will finally be output 
	#each "line" of the output refers to one position in the genome
	#it needs to be able to hold many attributes and they all need to be extendable into lists because some positions can be part of multiple features
	def __init__(self, position_info): #initialize the instance
		self.n = 0
		self.seqid = position_info['seqid']
		self.ref_pos = position_info['ref_pos']
		self.ref_base = position_info['ref_base']
		self.name = [position_info['name']]
		self.cds_pos = [position_info['cds_pos']]#these are currently zero-based
		self.strand = [position_info['strand']]
		self.frame = [position_info['frame']]
		self.codon = [position_info['codon']]
		self.aa = [position_info['aa']]
		self.degen = [position_info['degen']]
	def add(self, position_info):
		if self.seqid != position_info['seqid'] and self.ref_pos != position_info['ref_pos'] and ref_base != position_info['ref_base']:
			#some thing has gone wrong
			print("These aren't the same transcript, something has gone haywire!")
			sys.exit()
			return False
		else:
			self.n += 1
			self.name += [position_info['name']]
			self.cds_pos += [position_info['cds_pos']]
			self.strand += [position_info['strand']]
			self.frame += [position_info['frame']]
			self.codon += [position_info['codon']]
			self.aa += [position_info['aa']]
			self.degen += [position_info['degen']]
	def outstring(self):
		#the formatted string to out output to file/screen
		output_string = "\t".join(str(thing) for thing in [ \
		self.seqid, self.ref_pos, self.ref_base, \
		":".join(self.name), \
		":".join([str(i) for i in self.cds_pos]), \
		":".join(self.strand), \
		":".join([str(j) for j in self.frame]), \
		":".join(self.codon), \
		":".join(self.aa),\
		":".join([str(k) for k in self.degen]) \
		])
		return output_string
 

############################################################################################################################################
class Transcript:
	"""This is the core of the annotation
	The Transcript object holds all the information for a particular feature
	primarily, for a given transcript all its features are stored in a dictionary
	callled features. This dictionary has keys which refer to each feature type 
	and the value is another dictionary of all features of that type ordered by
	the start position of the feature"""
	#dict to look up AA's
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
	#dict to look up degeneracy
	degen_dict = {'ttt': '002', 'ttc': '002', 'tta': '202', 'ttg': '202', \
		'tct': '004', 'tcc': '004', 'tca': '004', 'tcg': '004', \
		'tat': '002', 'tac': '002', 'taa': '022', 'tag': '002', \
		'tgt': '002', 'tgc': '002', 'tga': '020', 'tgg': '000', \
		'ctt': '004', 'ctc': '004', 'cta': '204', 'ctg': '204', \
		'cct': '004', 'ccc': '004', 'cca': '004', 'ccg': '004', \
		'cat': '002', 'cac': '002', 'caa': '002', 'cag': '002', \
		'cgt': '004', 'cgc': '004', 'cga': '204', 'cgg': '204', \
		'att': '003', 'atc': '003', 'ata': '003', 'atg': '000', \
		'act': '004', 'acc': '004', 'aca': '004', 'acg': '004', \
		'aat': '002', 'aac': '002', 'aaa': '002', 'aag': '002', \
		'agt': '002', 'agc': '002', 'aga': '202', 'agg': '202', \
		'gtt': '004', 'gtc': '004', 'gta': '004', 'gtg': '004', \
		'gct': '004', 'gcc': '004', 'gca': '004', 'gcg': '004', \
		'gat': '002', 'gac': '002', 'gaa': '002', 'gag': '002', \
		'ggt': '004', 'ggc': '004', 'gga': '004', 'ggg': '004'}
	# init has replaced Dan's "new():"
	def __init__(self, l, index_label='pacid'): #####fragile
		"""This is a transcript object which deals with mRNA, exons, CDS etc
		The name of the transcript must be a unique identifier that is pulled from the attributes
		To initialize a transcript you iniitialize the class with a line from the GFF"""
		#l is a gff line object defined above
		#index_label is the thing in the attributes you want to use to uniquely identify transcripts 
		self.name = l.attributes[index_label]
		if 'Parent' in l.attributes:
			self.gene = l.attributes['Parent'] 
		self.fadb = ''  # this is Dan's fasta DB of the reference sequence - I will likely make a dict in the end script
		self.seqid = l.seqid #'seq_id' #this is the name of the fasta sequence in the reference
		self.start = l.start
		self.end = l.end
		self.strand = l.strand
		self.feats = {'mRNA': {}, \
					'exon': {}, \
					'CDS': {}, \
					'five_prime_UTR': {}, \
					'three_prime_UTR': {}, \
					'gene':{}}
		self.feats[l.type][l.start] = l
		self._cds = False
		self._cds_ref = False
		self._cds_length = False
		self._cds_map = False
		self._cds_frame = False
		self._valid_frame = False
		self._valid_cds = False
		self._aa = False
		self._cds_degen = False
		self._codons = False
	#this adds a new feature to a transcript
	def add(self, l, index_label='pacid'):
		if l.attributes[index_label] != self.name: 
			print('what the hell? Genes do not match')
			sys.exit()
		if l.start < self.start: self.start = l.start #this just sets the start and end of the feature as the first and last position of the feature
		if l.end > self.end: self.end = l.end
		if l.type in self.feats:
			self.feats[l.type][l.start] = l #the features of a transcript are kept in a list indexed by their start position
		else:
			self.feats[l.type]={l.start:l}
	#poorly named, this returns a list of all features of a given kind (CDS) for a given transcript, ordered by start position
	def sorted_feats(self, feature_type):
		#the feats dict of every transcript instance will have a number of dicts in it
		#each of these dicts will correspond to different types of features
		#each of these features will be a dict where the key is the start position of the feature and the value is the feature itself [GFF object]
		#this function will return all of the instances of any one kind of feat ordered by start position
		# I think this line is a problem because once this is called once it will never be re-run
		# this means that it will be run for exons and then not changed it you call it for CDS
		#if self._sorted_fragments: return self._sorted_fragments
		fragment_starts = list(self.feats[feature_type].keys())
		fragment_starts.sort()
		sorted_fragments = [self.feats[feature_type][i] for i in fragment_starts]
		#this is a list gff_line objects sorted by start position
		# These need to be removed to avoid this thing being stored
		# self._sorted_fragments = sorted_fragments
		# return self._sorted_fragments
		return sorted_fragments
	#returns CDS chunks all concatenated together
	def cds_ref(self, ref_fasta_dict):
		#this is the chunks of codings sequence pulled out of the reference, reversed if necessary but not complemented for frame
		if self._cds_ref: return self._cds_ref
		sorted_fragments = self.sorted_feats('CDS')
		cds_ref_seq = ''
		for frag in sorted_fragments:
			cds_ref_seq += ref_fasta_dict[frag.seqid].seq[frag.start-1:frag.end] #note GFF is 1-based and python is 0-based hence (start to end)  = [start-1:end]
		if self.strand == "-": cds_ref_seq = cds_ref_seq[::-1]
		self._cds_ref = cds_ref_seq
		return self._cds_ref
	#returns cds as a list reversed complementd and all buffed up
	def cds(self, ref_fasta_dict):
		if self._cds: return self._cds
		if self.strand == "-":
			coding_sequence = reverse_complement(self.cds_ref(ref_fasta_dict)[::-1])
		else:
			coding_sequence = self.cds_ref(ref_fasta_dict)
		self._cds = list(coding_sequence)
		return self._cds
	#calculates and stores the length of the CDS 
	def cds_length(self):
		if self._cds_length: return self._cds_length
		self._cds_length = len(self.cds_map())
		return self._cds_length
		#calculates and stores the length of the CDS 
	#list of genomic position of each position of the coding sequence (1-based)
	def cds_map(self):
		if self._cds_map: return self._cds_map
		sorted_fragments = self.sorted_feats('CDS')
		position_map = []
		for frag in sorted_fragments:
			position_map += range(frag.start,frag.end+1)
		if self.strand == "-":
			position_map = position_map[::-1]
		self._cds_map = position_map
		return self._cds_map
		#allows you to map positions in the CDS to positions in the genome reference
		#The numbers refer to the 1-based genome positions
		#I think its a list of numbers, each one a position in the genome, ordered by their appearance in the CDS (after RC etc)
		#eg	CDS	1	2	3	4	5	6	7	8	9	0	
		#	REF	100	101	102	103	210	211	212	213	350	351
		#	Exo	1	1	1	1	2	2	2	2	3	3
		#Note this will be reversed if necessary
#	def cds_exonno(self):
		# stores the exon number for each position in the CDS as an array
		# exon number is pulled from the GFF/GTF directly and not 'computed'
	def cds_frame(self):
		if self._cds_frame: return self._cds_frame
		#returns the frame of each position in the CDS as in 0,1,2
		#self._cds_frame = [0,1,2] * (self.cds_length()/3) # What does this do when the CDS isn't a multiple of 3!!! How do we know a CDS starts on position 1?
		self._cds_frame = ([0,1,2] * int((self.cds_length()+2)/3) )[:self.cds_length()] # DOES THIS FIX IT????
		return self._cds_frame
	def valid_frame(self):
		if self._valid_frame: return self._valid_frame
		# checks if the cds_frames match the frames encoded in GTF
		cds_position_counter = 0
		error = 0
		chunks = self.sorted_feats('CDS')
		if chunks[0].strand == "-": chunks = chunks[::-1]
		for chunk in chunks:
			GFF_phase = chunk.phase
			frame_trans = dict(zip('012', '021'))
			GFF_frame = "".join(frame_trans[i] for i in GFF_phase)
			# PYTHON 2 GFF_frame = GFF_phase.translate(maketrans('012', '021'))
			my_frame = str(self.cds_frame()[cds_position_counter])
			if my_frame != GFF_frame: 
				error+=1
			cds_position_counter += (chunk.end - chunk.start +1)
		if error > 0: 
			self._valid_frame = False
			return False
		else: 
			self._valid_frame = True
			return True 
		# each CDS chunk, the first position should have a phase
		# GFF_phase = feats['CDS'][start_position_1].phase #= 0, 1, or 2
		# my_frame =  
		# its cds_frame is the position that where
		# if phase = 0 frame should == 0
		# if phase = 1 frame should == 2
		# if phase = 2 frame should == 1
		# other than counting into the position map how do I know the frame of a given position
		# note that the frame col in a GTF file should be "phase" which is different!
		# see here for explanation: http://gmod.org/wiki/GFF
		# No. of bases removed from the start to find the beginning of the next full codon
		# eg. (ATG, phase 0, frame 1), (xATG, phase 1, frame 3), (xxATG, phase 2, frame 2)
	def aa(self, ref_fasta_dict):
		if self._aa: return self._aa
		codon_sequence = self.codons(ref_fasta_dict)
		aa_sequence = []
		for codon in codon_sequence:
			if codon in self.gen_code_dict:
				aa_sequence += [self.gen_code_dict[codon]]
			else:
				aa_sequence += ["."]
		self._aa = aa_sequence
		return self._aa
	#checks for internal starts and stops, and start and stop at ends
	def valid_cds(self, ref_fasta_dict):
		if self._valid_cds: return self._valid_cds
		aa = self.aa(ref_fasta_dict)
		if aa[-1] == 'X' and aa[0] == 'M' and 'X' not in aa[:-1]: 
			self._valid_cds = True
			return True
		else: 
			self._valid_cds = False
			return False
	def codons(self, ref_fasta_dict):
		#returns a list of the codons (ie length CDS/3)
		if self._codons:return self._codons
		coding_string = ''.join(self.cds(ref_fasta_dict)).lower()
		codon_sequence = [coding_string[i:i+3] for i in range(0, len(coding_string), 3)]
		self._codons = codon_sequence
		return self._codons
	def cds_degen(self, ref_fasta_dict):
		#return degeneracy for each position in the CDS as a list
		if self._cds_degen:return self._cds_degen
		codon_sequence = self.codons(ref_fasta_dict)
		degen_sequence = []
		for codon in codon_sequence: 
			if codon in self.degen_dict: degen_sequence += list(self.degen_dict[codon])
			else:degen_sequence += list("...")
		self._cds_degen = degen_sequence
		return self._cds_degen
	#def degen():made useless with my degen dict.
	def position_info(self, cds_pos, ref_fasta_dict):
		#takes a single CDS positions and returns a dict with some shit in it.
		pos_info = { \
		'name': self.name, \
		'seqid': self.seqid, \
		'cds_pos': cds_pos, \
		'ref_pos': self.cds_map()[cds_pos], \
		'ref_base': ref_fasta_dict[self.seqid][self.cds_map()[cds_pos]-1], \
		'strand': self.strand, \
		'frame': self.cds_frame()[cds_pos], \
		'codon': self.codons(ref_fasta_dict)[cds_pos/3], \
		'aa': self.aa(ref_fasta_dict)[cds_pos/3], \
		'degen': self.cds_degen(ref_fasta_dict)[cds_pos]}
		return pos_info
	def introns(self):
		"""this function will return a new feature type of 'introns'
		it will be inserted into the feats dictionary like other transcript features
		It is generated by finding the coordinates between exons
		Note exons include 5'UTR 3'UTR so using only CDS is inappropriate"""
		sorted_exons = self.sorted_feats('exon')
		if 'intron' in self.feats: print("introns already exist!")
		else:
			self.feats['intron'] = {}
			for i in range(len(sorted_exons)-1):
				upstream_exon =  sorted_exons[i]
				downstream_exon = sorted_exons[i+1]
				intron_line = "\t".join([str(j) for j in [\
				self.seqid, 'inferred', 'intron', upstream_exon.end+1, downstream_exon.start-1, \
				".", self.strand, ".", "ness_ID=%s;intron_number=%i" %(self.name, i+1)]])
				intron_l = GFF_line.GFF_line(intron_line)
				self.feats['intron'][intron_l.start] = intron_l
		return self.feats['intron']


def hash_gff(gff_file, index_label = 'pacid',info_delimiter=";", info_field_delimiter = '=', quiet=False):
	"""This function will take a GFF and make a dictionary of transcript objects.
	It requires an index label that uniquely identifies a transcript
	"""
	gff=open(gff_file, 'r')
	transcripts = {}
	for line in gff:
		if line[0] == "#": pass
		else:
			g = GFF_line.GFF_line(line.strip(), info_delimiter, info_field_delimiter)
			if index_label in g.attributes:
				if g.attributes[index_label] not in transcripts:
					transcripts[g.attributes[index_label]] = Transcript(g, index_label=index_label)
				else:
					transcripts[g.attributes[index_label]].add(g, index_label=index_label)
			elif not quiet:
				print("this feature does not have an appropriate index")
				print(line)
	gff.close()
	return transcripts
