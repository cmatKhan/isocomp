import sys
from importlib.metadata import version 
from .create_windows import main as create_windows
from .rename_fasta import main as rename_fasta


#--------------------------------------------------------
# parser call_variants
#--------------------------------------------------------
# positional/optional arguments
posList = ['input', 'output']
optList = ['N'] #
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-30)'
description = '\nDescription: The program calls SNPs and InDels from the unique isoform file'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput file.  e.g. /in/unique_isoforms.tsv')
parser.add_argument(posList[1], type=str, help='string\toutput file.  e.g. /out/unique_isoforms_variants.tsv')
# optional arguments
# NOTE: currently called N. renamed to max_indels
parser.add_argument('-n', '--'+optList[0], type=int, metavar='', default=10, help='float\tmaximum bases for indel combination (inclusive), default 10')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()
for key, value in argsDict.items():
    if key in posList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    if key in optList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Optional Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    vars()[key] = value # assign values of arguments into shorthand global variables
sys.stdout.write('\n'); sys.stdout.flush()

for key, value in argsDict.items():
    if key in posList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    if key in optList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Optional Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    vars()[key] = value # assign values of arguments into shorthand global variables
sys.stdout.write('\n'); sys.stdout.flush()

#--------------------------------------------------------
# parser coverage stats
#--------------------------------------------------------
# positional/optional arguments
posList = ['bam', 'bed', 'samtools', 'output']
optList = ['mapQCut']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-11)'
description = '\nDescription: The program generate bam coverage statistics for given bed regions'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput bam file.  e.g. /in/sample.bam')
parser.add_argument(posList[1], type=str, help='string\tinput bed file.  e.g. /in/sample.bed')
parser.add_argument(posList[2], type=str, help='string\tsamtools.  e.g. /path/to/samtools')
parser.add_argument(posList[3], type=str, help='string\toutput regional coverage statistics.  e.g. /out/sample.bed')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=int, metavar='', default=20, help='integer\tminimum mapQ for counted reads, default 20')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()


def main()->None:
	"""Entry point to isocomp"""

	isocomp_version = version('isocomp')

	version_stmnt = f'Isocomp version: {isocomp_version}'  
	# Stuff to print before offering a list of tools
	preamble = [
		version_stmnt,
		'Description: A set of tools to discover novel isoforms',
		'Available tools are:',
	]

    # list of available tools
	# preface each with \t
	tool_list = [
		'\tcreate_windows',
		'\trename_fasta'
	]

	# Anything to print after the list of tools
	epigraph = [
		'For usage instructions, enter isocomp <tool_name> --help'
	]

	help_str = "\n".join(
		["\n".join(preamble),"\n".join(tool_list),"\n".join(epigraph)])
	
	try:
		tool = sys.argv[1]
	except IndexError:
		tool = 'default'

	if tool == '--version': 
		print(version_stmnt)
	elif tool == '--help': 
		print(help_str)
	elif tool == 'create_windows': 
		create_windows(sys.argv[2:])
	elif tool == 'rename_fasta':
		rename_fasta(sys.argv[2:])
	else:
		print(help_str)

if __name__ == "__main__":
    sys.exit(main())
