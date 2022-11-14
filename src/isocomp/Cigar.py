# pylint:disable=W0622,C0103
# standard library
import re
import logging
# outside dependencies

__all__ = ['Cigar']

logging.getLogger(__name__).addHandler(logging.NullHandler())

# TODO document class
class Cigar:
	"""_summary_"""

	def __init__(self):
		self._string = ''  # cigar string
		self._ref = ''  # ref sequence
		self._query = ''  # query sequence
		self._ref_id = ''  # ref id
		self._query_id = ''  # query_id
		self._parsed = []  # parsed cigar
		self._snps = []  # called snps
		self._indels = []  # called indels

	# TODO document attributes.
	# see https://docs.python.org/3.9/library/functions.html#property
	@property
	def string(self):
		"""_summary_"""
		return self._string

	@string.setter
	def string(self, new_string):
		self._string = new_string

	@property
	def ref(self):
		"""_summary_"""
		return self._ref

	@ref.setter
	def ref(self, new_ref):
		self._ref = new_ref

	@property
	def query(self):
		"""_summary_"""
		return self._query

	@query.setter
	def query(self, new_query):
		self._query = new_query

	@property
	def ref_id(self):
		"""_summary_"""
		return self._ref_id

	@ref_id.setter
	def ref_id(self, new_ref_id):
		self._ref_id = new_ref_id

	@property
	def query_id(self):
		"""_summary_"""
		return self._query_id

	@query_id.setter
	def query_id(self, new_query_id):
		self._query_id = new_query_id

	@property
	def parsed(self):
		"""_summary_"""
		return self._parsed

	@parsed.setter
	def parsed(self, new_parsed):
		self._parsed = new_parsed

	@property
	def snps(self):
		"""_summary_"""
		return self._snps

	@snps.setter
	def snps(self, new_snps):
		self._snps = new_snps

	@property
	def indels(self):
		"""_summary_"""
		return self._indels

	@indels.setter
	def indels(self, new_indels):
		self._indels = new_indels

	# function to call snps and indels
	def call(self):

		self.parsed = re.findall(r'(\d+)(\w)', self.string.replace('=', 'M'))
		self.parsed = [[int(c[0]), c[1]] for c in self.parsed]

		ref_pos, query_pos = 1, 1

		for c in self.parsed:

			# matches
			if c[1] == 'M':
				ref_pos += c[0]
				query_pos += c[0]
			# mismatches
			elif c[1] == 'X':
				self.snps.append(
					['snp',
					 ref_pos,
					 ref_pos+c[0],
					 query_pos,
					 query_pos+c[0],
					 self.ref[ref_pos-1:ref_pos+c[0]-1],
					 self.query[query_pos-1:query_pos+c[0]-1]])
				ref_pos += c[0]
				query_pos += c[0]
			# deletions in query
			elif c[1] == 'D':
				self.indels.append(
					['del',
					 ref_pos,
					 ref_pos+c[0],
					 query_pos,
					 query_pos,
					 self.ref[ref_pos-1:ref_pos+c[0]-1], '.'])
				ref_pos += c[0]
			# insertions in query
			elif c[1] == 'I':
				self.indels.append(
					['ins',
					 ref_pos,
					 ref_pos,
					 query_pos,
					 query_pos+c[0],
					 '.', self.query[query_pos-1:query_pos+c[0]-1]])
				query_pos += c[0]
			else:
				raise ValueError(f"{c[0]} is not a recognized cigar notation")
    
	# TODO document input -- check input type
	def refine(self, cut:int) -> None:
		"""merge close indels and remove condensed snps

		Args:
			cut (int): _description_
		"""
		update_snps, update_indels, last = [], [], []

		# combine close ins or del
		for current in self.indels:

			# handle the starting case
			if last == []:
				update_indels.append(current)
				last = update_indels[-1]
				continue

			# retrieve the last indel
			last = update_indels[-1]
			type_l, ref_start_l, ref_end_l, \
				query_start_l, query_end_l, ref_l, alt_l = last
			type_c, ref_start_c, ref_end_c, \
				query_start_c, query_end_c, ref_c, alt_c = current

			# current indel to be merged with the last one
			if (type_l == type_c or type_l == 'delins') and \
					ref_start_c-ref_end_l <= cut and \
					query_start_c-query_end_l <= cut:
				update_indels[-1] = \
					['delins', 
					 ref_start_l, 
					 ref_end_c, 
					 query_start_l, 
					 query_end_c,
					 self.ref[ref_start_l-1:ref_end_c-1], 
					 self.ref[query_start_l-1:query_end_c-1]]
			# current indel not to be merged
			else:
				update_indels.append(current)

		# check snps
		for snp in self.snps:
			for indel in update_indels:

				# snps contained in a merged delins
				typeS, ref_start_s, ref_end_s, \
					query_start_s, query_end_s, ref_s, alt_s = snp
				typeI, ref_startI, ref_end_i, \
					query_start_i, query_end_i, ref_i, alt_i = indel

				if ref_start_s >= ref_startI and \
					ref_end_s <= ref_end_i and \
						query_start_s >= query_start_i and \
							query_end_s <= query_end_i:
					update_snps.append(snp)
					break

		self.snps = [s for s in self.snps if s not in update_snps]
		self.indels = update_indels
