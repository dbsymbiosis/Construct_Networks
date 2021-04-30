#!/usr/bin/env python2
DESCRIPTION = '''
Takes the node annotation file (i.e. KEGG_Pathway_Networks.nodes_withSeqIDs.txt.gz) after sequence and compound IDS
have been added, and adds weather or not a sequence or gene annotated to that node has been identified as differentially
expressed or accumulated. 

Will also added an extra column detailing info about which condition/s each sequence or compound was found to be 
differentially expressed or accumulated under. This will be helpful for downstream analysis

NOTE:
	- 
'''
import sys
import os
import argparse
import logging
import gzip

## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-n', '--nodes', metavar='KEGG_Pathway_Networks.nodes_withSeqIds.txt.gz',
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Nodes to annotate with info about diff. expression or accumulation (default: stdin)'
	)
	parser.add_argument('--diff_expr', metavar='diff_expr_genes.txt',
		required=False, default=None, type=lambda x: File(x, 'r'),
		help='Input gene_id[<tab>cond_info] mapping [gzip] file'
	)
	parser.add_argument('--diff_accum', metavar='diff_accum_metabolites.txt',
		required=False, default=None, type=lambda x: File(x, 'r'),
		help='Input compound_id[<tab>cond_info] mapping [gzip] file'
	)
	parser.add_argument('-o', '--out', metavar='KEGG_Pathway_Networks.nodes_withSeqIds_DiffExprAccum.txt.gz',
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'),
		help='Output [gzip] file with nodes annotated with info about diff. expression or accumulation (default: stdout)'
	)
	parser.add_argument('--debug', 
		required=False, action='store_true', 
		help='Print DEBUG info (default: %(default)s)'
	)
	args = parser.parse_args()
	
	## Set up basic debugger
	logFormat = "[%(levelname)s]: %(message)s"
	logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	
	logging.debug('%s', args) ## DEBUG
	
	diff_expr = {}
	diff_accum = {}
	if args.diff_expr is not None:
		with args.diff_expr as diff_expr_fh:
			diff_expr = load_id2gene(diff_expr_fh)
	if args.diff_accum is not None:
		with args.diff_accum as diff_accum_fh:
			diff_accum = load_id2gene(diff_accum_fh)
	
	with args.nodes as infile, args.out as outfile:
		add_diffExprAccum_to_Nodes(infile, outfile, diff_expr, diff_accum)



def add_diffExprAccum_to_Nodes(infile, outfile, diff_expr, diff_accum, col_delim='\t', id_delim=';'):
	'''
	Read input node annotation file and adds info about weather and which sequecnes or compounds are differentially 
	expressed or accumulated. 
	
	NOTE:
		- Assumes first line is header. 
		- Used header fow to find correct columns for analysis.
	'''
	header_line = infile.readline().strip('\n')
	headers = header_line.split(col_delim)
	
	gene_compound_ids_index = get_header_index(headers, "gene-compound_ids")
	type_index = get_header_index(headers, "type")
	
	outfile.write(header_line + col_delim + col_delim.join(['diff_expr-accum', 'diff_expr-accum_info_1', 'diff_expr-accum_info_2']) + '\n')
	for line in infile:
		line = line.strip('\n')
		if not line or line.startswith('#'):
			continue
		
		## Split and get values using index
		line_split = line.split(col_delim)
		gene_compound_ids_value = get_value_using_index(line_split, gene_compound_ids_index).split(id_delim)
		gene_compound_ids_value = [x for x in gene_compound_ids_value if x != "missing"] # Remove "missing" place holder values
		type_value = get_value_using_index(line_split, type_index)
		
		## Get seq/compund ids associated with node. Check which type of node we are looking and and search the right lists depending. 
		has_annots = False
		info_1 = []
		info_2 = []
		if type_value == "reaction":
			for i in gene_compound_ids_value:
				has_annots = True
				if i in diff_expr.keys():
					info_1.append(i+":"+":".join(diff_expr[i]))
					info_2.extend(diff_expr[i])
		elif type_value == "compound":
			for i in gene_compound_ids_value:
				has_annots = True
				if i in diff_accum.keys():
					info_1.append(i+":"+":".join(diff_accum[i]))
					info_2.extend(diff_accum[i])
		
		## If we didnt find any diff. expr or accum features associated with this node add missing to the two columns.
		if has_annots:
			if len(info_1) == 0:
				diff = 'No'
				info_1 = ['-']
				info_2 = ['-']
			else:
				diff = 'Yes'
		else:
			diff = 'missing'
			info_1 = ['-']
			info_2 = ['-']
		
		## Write new row to output file
		outfile.write(col_delim.join(line_split) + col_delim + col_delim.join([diff, '---'.join(info_1), ';'.join(set(info_2))]) + '\n')


def get_value_using_index(row_list, index):
	'''
	Takes a row split into columns and returns the value using the index provided.
	Will return an error is the index is missing from the row_list.
	'''
	try:
		return row_list[index]
	except IndexError:
		logging.error('Index "%s" is missing from row:\n%s', index, row_list) ## ERROR
		sys.exit(1)
	
	

def get_header_index(headers, header2find):
	'''
	Takes a list of headers and returns the index of the desired header.
	Will return an error if the header is missing from the list of headers. 
	'''
	count = headers.count(header2find)
	if count == 0:
		logging.error('Header "%s" is missing from header row:\n%s', header2find, headers) ## ERROR
		logging.error('Can\'t continue with this header missing - Stopping!') ## ERROR
		sys.exit(1)
	elif count == 1:
		idx = headers.index(header2find)
		logging.debug('Index for header "%s" is %s', header2find, idx) ## DEBUG
		return idx
	else:
		logging.error('Header "%s" is present %s times in the header row:\n%s', header2find, count, headers) ## ERROR
		logging.error('Unsure which one to pick - Stopping!') ## ERROR
		sys.exit(1)



def load_id2gene(fh, delim='\t'):
	'''
	Takes a input 2 column mapping file (key_id<tab>value) and returns a dictionary 
	of the mappings. 
	
	NOTE:
		- Does not remove duplicate values, assumes multiple is significant and that they can be filtered out later on.
	'''
	id2ids = {}
	for line in fh:
		line  = line.strip()
		if not line or line.startswith('#'):
			continue
		key, value = line.split(delim)
		if key not in id2ids.keys():
			id2ids[key] = []
		id2ids[key].append(value)
	return id2ids


class File(object):
	'''
	Context Manager class for opening stdin/stdout/normal/gzip files.

	 - Will check that file exists if mode='r'
	 - Will open using either normal open() or gzip.open() if *.gz extension detected.
	 - Designed to be handled by a 'with' statement (other wise __enter__() method wont 
	    be run and the file handle wont be returned)
	
	NOTE:
		- Can't use .close() directly on this class unless you uncomment the close() method
		- Can't use this class with a 'for' loop unless you uncomment the __iter__() method
			- In this case you should also uncomment the close() method as a 'for'
			   loop does not automatically cloase files, so you will have to do this 
			   manually.
		- __iter__() and close() are commented out by default as it is better to use a 'with' 
		   statement instead as it will automatically close files when finished/an exception 
		   occures. 
		- Without __iter__() and close() this object will return an error when directly closed 
		   or you attempt to use it with a 'for' loop. This is to force the use of a 'with' 
		   statement instead. 
	
	Code based off of context manager tutorial from: https://book.pythontips.com/en/latest/context_managers.html
	'''
 	def __init__(self, file_name, mode):
		## Upon initializing class open file (using gzip if needed)
		self.file_name = file_name
		self.mode = mode
		
		## Check file exists if mode='r'
		if not os.path.exists(self.file_name) and mode == 'r':
			raise argparse.ArgumentTypeError("The file %s does not exist!" % self.file_name)
	
		## Open with gzip if it has the *.gz extension, else open normally (including stdin)
		try:
			if self.file_name.endswith(".gz"):
				#print "Opening gzip compressed file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = gzip.open(self.file_name, self.mode+'b')
			else:
				#print "Opening normal file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = open(self.file_name, self.mode)
		except IOError as e:
			raise argparse.ArgumentTypeError('%s' % e)
	def __enter__(self):
		## Run When 'with' statement uses this class.
		#print "__enter__: %s" % (self.file_name) ## DEBUG
		return self.file_obj
	def __exit__(self, type, value, traceback):
		## Run when 'with' statement is done with object. Either because file has been exhausted, we are done writing, or an error has been encountered.
		#print "__exit__: %s" % (self.file_name) ## DEBUG
		self.file_obj.close()
#	def __iter__(self):
#		## iter method need for class to work with 'for' loops
#		#print "__iter__: %s" % (self.file_name) ## DEBUG
#		return self.file_obj
#	def close(self):
#		## method to call .close() directly on object.
#		#print "close: %s" % (self.file_name) ## DEBUG
#		self.file_obj.close()


if __name__ == '__main__':
	main()
