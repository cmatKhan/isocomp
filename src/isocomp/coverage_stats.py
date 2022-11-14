import os
import logging
from argparse import Namespace

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['coverage_stats']

# TODO better documentation on summary, args and returns 
def configure_bedfile(in_bed:str) -> dict:
	"""config input bed file. Not exported to isocomp API.

	Args:
		in_bed (str): _description_

	Raises:
		FileNotFoundError: if in_bed does not exist

	Returns:
		dict: _description_
	"""
	# note that the in_bed is checked in coverage_stats
 
	bed_dict = {}
	with open(in_bed) as f:

		for line in f:
			chr,start,end = line.strip().split()[:3]	
			#if chr in bed_dict: bed_dict[chr].append([start,end])
			#if not chr in bed_dict: bed_dict[chr] = [[start,end]]
			bed_dict.setdefault(chr,[]).append([start,end])

	return(bed_dict)

def coverage_stats(args:Namespace) -> None:

	# Check input paths
	logging.debug(args)

	input_path_list = [args.bed]
	for input_path in input_path_list:
		if not os.path.exists(input_path):
			error_msg=f"Input file DNE: {input_path}"
			logging.error(error_msg)
			raise FileNotFoundError(error_msg)

	# coverage statistics over target regions

	bed_dict = configure_bedfile(args.bed)

	out = open(args.output, 'w')
	#out.write( '\t'.join( ['chr','start','end','length','depthTotal','depthAve'] ) + '\n' ) # header for args.output

	for chr,beds in bed_dict.items():

		for bed in beds:

			start,end = bed
			logging.info('handling bed %s:%s-%s' %(chr,start,end))

			# obtain depth stats by samtools depth
			# TODO this needs to be done via pysam. pysam handles the samtools
			# dependency
			os.system('%s depth %s -r %s:%s-%s -Q %s > %s.tmp' %(samtools,args.bam,chr,start,end,mapQCut,args.output))
			depthDict = {}
			with open(args.output+'.tmp') as f:
				for line in f:
					depthDict[ int(line.strip().split('\t')[1]) ] = int(line.strip().split('\t')[2])

			# statistics
			for head in range(int(start),int(end),100):
				tail = head + 99
				if tail > int(end): tail = int(end)
				segTotal = sum( [ depthDict[j] for j in range(head,tail+1) if j in depthDict ] )
				segLength = tail - head + 1
				segAve = segTotal / segLength

				out.write( '\t'.join( map(str,[chr,head,tail,format(segAve,'.2f')]) ) + '\n' )
			# TODO use a context manager and temp directory via python stdlib 
			# temp rather than rm'ing directly
			# https://docs.python.org/3/library/tempfile.html#tempfile-examples
			os.system('rm %s.tmp' %(args.output))

	out.close()

# move to logger in __main__
# sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

