#!/usr/bin/env python2
DESCRIPTION = '''
Takes the node annotation file (i.e. KEGG_Pathway_Networks.nodes.txt.gz) and added the 
sequence IDs that have been annotated to each reaction and/or KEGG ortholog + the compound IDs
that have been annotated to each reaction.

NOTE:
	- Returns a list of unique ids per node (i.e. will remove duplicates if they exist in the annotations file)
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
	parser.add_argument('-n', '--nodes', metavar='KEGG_Pathway_Networks.nodes.txt.gz',
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Nodes to annotate with sequence IDs (default: stdin)'
	)
	parser.add_argument('--reaction2gene', metavar='reaction2gene.txt',
		required=False, default=None, type=lambda x: File(x, 'r'),
		help='Input reaction_id<tab>gene_id mapping [gzip] file'
	)
	parser.add_argument('--ortholog2gene', metavar='reaction2gene.txt',
		required=False, default=None, type=lambda x: File(x, 'r'),
		help='Input ortholog_id<tab>gene_id mapping [gzip] file'
	)
	parser.add_argument('--compound2gene', metavar='compound2gene.txt',
		required=False, default=None, type=lambda x: File(x, 'r'),
		help='Input compound_id<tab>gene_id mapping [gzip] file'
	)
	parser.add_argument('-o', '--out', metavar='KEGG_Pathway_Networks.nodes_withSeqIds.txt.gz',
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'),
		help='Output [gzip] file with nodes annotated with sequence or compound IDs (default: stdout)'
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
	
	reaction2gene = {}
	ortholog2gene = {}
	compound2gene = {}
	if args.reaction2gene is not None:
		with args.reaction2gene as reaction2gene_fh:
			reaction2gene = load_id2gene(reaction2gene_fh)
	if args.ortholog2gene is not None:
		with args.ortholog2gene as ortholog2gene_fh:
			ortholog2gene = load_id2gene(ortholog2gene_fh)
	if args.compound2gene is not None:
		with args.compound2gene as compound2gene_fh:
			compound2gene = load_id2gene(compound2gene_fh)
	
	with args.nodes as infile, args.out as outfile:
		add_seq_annots_to_Nodes(infile, outfile, reaction2gene, ortholog2gene, compound2gene)



def add_seq_annots_to_Nodes(infile, outfile, reaction2gene, ortholog2gene, compound2gene, col_delim='\t', id_delim=';'):
	'''
	Read input node annotation file and add reaction, ortholog, and compound ID seq and compound ids to them. 
	
	NOTE:
		- Assumes first line is header. 
		- Used header fow to find correct columns for analysis.
	'''
	header_line = infile.readline().strip('\n')
	headers = header_line.split(col_delim)
	kegg_id_index = get_header_index(headers, "kegg_id")
	info_index = get_header_index(headers, "info")
	type_index = get_header_index(headers, "type")
	
	outfile.write(header_line + col_delim + 'gene-compound' + col_delim + 'gene-compound_ids\n')
	for line in infile:
		line = line.strip('\n')
		if not line or line.startswith('#'):
			continue
		
		## Split and get values using index
		line_split = line.split(col_delim)
		kegg_id_value = get_value_using_index(line_split, kegg_id_index).split(id_delim)
		info_value = get_value_using_index(line_split, info_index).split(id_delim)
		type_value = get_value_using_index(line_split, type_index)
		
		## Get seq/compund ids associated with node. Check which type of node we are looking and and search the right lists depending. 
		ids_found = []
		if type_value == "reaction":
			for i in kegg_id_value:
				if i in reaction2gene.keys():
					ids_found.extend(reaction2gene[i])
			for i in info_value:
				if i in ortholog2gene.keys():
					ids_found.extend(ortholog2gene[i])
		elif type_value == "compound":
			for i in kegg_id_value:
				if i in compound2gene.keys():
					ids_found.extend(compound2gene[i])
		elif type_value == "ortholog":
			for i in kegg_id_value:
				if i in ortholog2gene.keys():
					ids_found.extend(ortholog2gene[i])
		elif type_value == "gene":
			for i in info_value:
				if i in ortholog2gene.keys():
					ids_found.extend(ortholog2gene[i])
		
		## If we didnt find any annotations add missing to the two columns. 
		if len(ids_found) == 0:
			ids_found = ['missing']
			ids = 'missing'
		else:
			ids = 'present'
		
		## Write new row to output file
		outfile.write(col_delim.join(line_split) + col_delim + ids + col_delim + id_delim.join(set(ids_found)) + '\n')


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
