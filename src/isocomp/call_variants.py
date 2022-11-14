"""
call_variants unique_isoforms_and_aln_stats_filter0.99.tsv unique_isoforms_and_aln_stats_filter0.99.call.tsv --N 10



"""
import os
import logging
from argparse import Namespace

from .Cigar import Cigar

logging.getLogger(__name__).addHandler(logging.NullHandler())

# TODO Docstring summary. format the output so that it displays as plaintext
# in mkdocs. better args description -- probably best to put example of the 
# Namespace to pass. Usage example is in the parse_args in __main__
def call_variants(args:Namespace) -> None:
	"""_summary_
	
	## Output
	ref_id: reference sequence ID
	query_id: query sequence ID
	type: type of variant (snp/ins/del/delins)
	ref: variant sequence in reference
	alt: altered sequence
	ref_start: starting position of the variant in reference

	#ref_id	query_id	type	ref	alt	ref_pos
	HG002_PB.6.2	HG004_PB.6.2	delins	G	GCCAAGGTCC	4
	HG002_PB.6.2	HG004_PB.6.2	ins	.	AAAGATATTTC	1291
	HG002_PB.6.2	HG005_PB.7566.3	snp	A	G	7

	Args:
		args (argparse.Namespace): Input from cmd line parser. see __main__
	"""
	
	# Check input paths
	logging.debug(args)

	input_path_list = [args.input]
	for input_path in input_path_list:
		if not os.path.exists(input_path):
			error_msg=f"Input file DNE: {input_path}"
			logging.error(error_msg)
			raise FileNotFoundError(error_msg)

	out = open(args.output, 'w')
	out.write('\t'.join(['#ref_id','query_id','type','ref','alt','ref_pos']) +'\n')

	with open(args.input) as f:

		used = [('','')]*100

		for index,line in enumerate(f):

			if index == 0: continue

			# parse input, skip entries with total==1
			try:
				chr, start, end, total, query, \
					sample_from, sample_to, map_start, \
						querySeq, alns = line.strip().split('\t')
				alns = alns.split(',')
			except: #pylint:disable=W0702
				continue

			# skip duplicated entries
			if (sample_from, query) in used: continue
			used = used[1:]
			used.append((sample_from, query))

			# obtain all candidate isoform alignments
			cigars = [ a.split('__')[-1] for a in alns ]
			ref_seqs = [ a.split('__')[-2] for a in alns ]
			ref_ids = [ a.split('__')[-3] for a in alns ]

			# parse cigar and call variants
			for i,c in enumerate(cigars):

				target = Cigar()
				target.string = c
				target.ref = ref_seqs[i]
				target.query = querySeq
				target.ref_id = ref_ids[i]
				target.query_id = sample_from+'_'+query

				# call variants
				target.call()
				# refine indels
				target.refine(args.max_indels)

				for v in target.snps + target.indels:
					out.write('\t'.join(
						map(str, 
						[target.ref_id,
						target.query_id,
						v[0],v[-2],v[-1],v[1]])) +'\n')




# move this to logger in __main__ sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

