from pysam import Tabixfile


#########################################################################################################################################
def reverse_complement(sequence):
	t = dict(zip('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB') )
	return "".join([t[b] for b in str(sequence)[::-1]])
############################################################################################################################################
class annotation_line(object):
	"""
	This is a simple class that allows easy access to the annotation table lines
	It is different from annotation_table_position in that it doesn't fetch the line from the table - it takes a line from the table and parses it up
	"""
	def __init__(self, line):
		self.line=line
		self.chromosome, \
		self.position, \
		self.reference_base, \
		self.genic, \
		self.exonic, \
		self.intronic, \
		self.intergenic, \
		self.utr5, \
		self.utr3, \
		self.fold0, \
		self.fold4, \
		self.fold2, \
		self.fold3, \
		self.CDS, \
		self.mRNA, \
		self.rRNA, \
		self.tRNA, \
		self.feature_names, \
		self.feature_types, \
		self.feature_ID, \
		self.cds_position, \
		self.strand, \
		self.frame, \
		self.codon, \
		self.aa, \
		self.degen, \
		self.FPKM, \
		self.rho, \
		self.FAIRE, \
		self.recombination, \
		self.mutability, \
		self.quebec_alleles = self.line.split('\t')
		self.position = int(self.position)


#####
class annotation_table_position(object):
	"""
	uses the genome annotation table I generated to retrieve a bunch of descriptive
	features of the genome
	
	I added a header line which could be used to define the features that way I don't
	have to modify this everytime I add a new column
	"""
	def __init__(self, chromosome, position, annotation_table_file ):
		annotation_table = Tabixfile(annotation_table_file)
		self.line = annotation_table.fetch(reference=chromosome, start=position-1, end=position).next()
		self.chromosome, \
		self.position, \
		self.reference_base, \
		self.genic, \
		self.exonic, \
		self.intronic, \
		self.intergenic, \
		self.utr5, \
		self.utr3, \
		self.fold0, \
		self.fold4, \
		self.fold2, \
		self.fold3, \
		self.CDS, \
		self.mRNA, \
		self.rRNA, \
		self.tRNA, \
		self.feature_names, \
		self.feature_types, \
		self.feature_ID, \
		self.cds_position, \
		self.strand, \
		self.frame, \
		self.codon, \
		self.aa, \
		self.degen, \
		self.FPKM, \
		self.rho, \
		self.FAIRE, \
		self.recombination, \
		self.mutability, \
		self.quebec_alleles = self.line.split('\t')
		self.position = int(self.position)
		annotation_table.close()
	
	def output_line(self):
		o_string = "\t".join( 
		[str(i) for i in \
		[self.chromosome, \
		self.position, \
		self.reference_base, \
		self.genic, \
		self.exonic, \
		self.intronic, \
		self.intergenic, \
		self.utr5, \
		self.utr3, \
		self.fold0, \
		self.fold4, \
		self.fold2, \
		self.fold3, \
		self.CDS , \
		self.mRNA, \
		self.rRNA, \
		self.tRNA, \
		self.feature_names, \
		self.feature_types, \
		self.feature_ID, \
		self.cds_position, \
		self.strand, \
		self.frame, \
		self.codon, \
		self.aa, \
		self.degen, \
		self.FPKM, \
		self.rho, \
		self.FAIRE, \
		self.recombination, \
		self.mutability, \
		self.quebec_alleles]])
		return o_string

############################################################################################################################################