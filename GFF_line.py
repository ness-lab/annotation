
import sys, re
from collections import OrderedDict

############################################################################################################################################
class GFF_line:
	"""This class basically parses a GFF line and allows you to interact with different components that I have deemed interesting
	Most components are simple strings or intgers.
	The attributes field which is a ;-separated list is returned as a dictionary """
	def __init__(self, l, info_delimiter=";", info_field_delimiter = '='):
		self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attribs = l.split('\t')
		self.attributes = self.attribute_dict(self.attribs, info_delimiter, info_field_delimiter)
		self.start = int(self.start)
		self.end = int(self.end)
		self.line = l
	def attribute_dict(self, attributes, info_delimiter=";", info_field_delimiter = '='): ###this is fragile
		d = OrderedDict()
		attributes = attributes.strip(info_delimiter)
		for i in [x.strip() for x in re.split(info_delimiter, attributes)]:
			if len(i)>0:
				field = re.split(info_field_delimiter, i)[0]
				d[field] = re.split(info_field_delimiter, i)[1].strip('"')
		return d
	def retrieve_sequence(self, ref_dict, reverse_complement=False):
		seq = str(ref_dict[self.seqid].seq[self.start-1: self.end])
		if reverse_complement:
			seq = reverse_complement(seq)
		return seq



def reverse_complement(sequence):
	#python 2 return str(sequence)[::-1].translate(maketrans('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))
	tr = dict(zip('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))
	return "".join([tr[i] for i in sequence])[::-1]