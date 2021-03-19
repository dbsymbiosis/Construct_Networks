#!/usr/bin/env python2
DESCRIPTION = '''
Filter MAGI gene results found in file magi_gene_results.csv
Filter by compound_score >= 1, reciprocal_score = 2, e_score_r2g > 5, e_score_g2r > 5

## Output (no header):
gene_id [tab] database_id_g2r

gene_id: ID of gene that has been annotated
		-  A single gene_id can have multiple annotations (multiple lines in outfile file)
database_id_g2r: ID of reaction annotation (only one ID per line)
		- Either RHEA or MetaCyc ID
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
	parser.add_argument('-i', '--input', metavar='magi_gene_results.csv', 
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Input [gzip] MAGI gene annotation (default: stdin)'
	)
	parser.add_argument('-o', '--out', metavar='magi_gene_results.filtered.txt', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] filtered annotated genes (default: stdout)'
	)
	parser.add_argument('--compound_score', 
		required=False, default=1, type=int, 
		help='Keep genes with compound_score >= X (default: %(default)s)'
	)
	parser.add_argument('--reciprocal_score',
                required=False, default=2, type=int,
                help='Keep genes with reciprocal_score == X (default: %(default)s)'
        )
	parser.add_argument('--e_score_r2g',
                required=False, default=5, type=int,
                help='Keep genes with e_score_r2g > X (default: %(default)s)'
        )
	parser.add_argument('--e_score_g2r',
                required=False, default=5, type=int,
                help='Keep genes with e_score_g2r > X (default: %(default)s)'
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
	
	
	with args.input as infile, args.out as outfile:
		filter_magi_gene_results(infile, outfile, args.compound_score, args.reciprocal_score, args.e_score_r2g, args.e_score_g2r)
	
	


def filter_magi_gene_results(infile, outfile, threshold_compound_score, threshold_reciprocal_score, threshold_e_score_r2g, threshold_e_score_g2r):
	'''
	Filter: compound_score >= X, reciprocal_score = X, e_score_r2g > X, e_score_g2r > X
	
	## Header/columns: magi_gene_results.csv
	# 1 MAGI_score
	# 2 gene_id
	# 3 original_compound
	# 4 neighbor
	# 5 note
	# 6 compound_score
	# 7 level
	# 8 homology_score
	# 9 reciprocal_score
	# 10 reaction_connection
	# 11 e_score_r2g
	# 12 database_id_r2g
	# 13 e_score_g2r
	# 14 database_id_g2r
	'''
	col_names = ["gene_id", "compound_score", "reciprocal_score", "e_score_r2g", "e_score_g2r", "database_id_g2r"]
	for gene_id, compound_score, reciprocal_score, e_score_r2g, e_score_g2r, database_id_g2r in iter_columns_by_name(infile, col_names):
		try:
			reciprocal_score = float(reciprocal_score)
			
			## If gene <-> compound annotation is not reciprocal then these values could be empty.
			## Set missing scores to -1 so we have something to evaluate later on
			## Since by default we only want reciprocal_score == 2 the missing values we use will have no actual affect downstream as none of these rows will pass by default. 
			missing = -1
			if compound_score != "":
				compound_score = float(compound_score)
			else:
				compound_score = missing
			
			if e_score_r2g != "":
				e_score_r2g = float(e_score_r2g)
			else:
				e_score_r2g = missing
			
			if e_score_g2r != "":
				e_score_g2r = float(e_score_g2r)
			else:
				e_score_g2r = missing
		except ValueError as e:
			logging.error('gene_id:"%s" compound_score:"%s" reciprocal_score:"%s" e_score_r2g:"%s" e_score_g2r:"%s"', gene_id, compound_score, reciprocal_score, e_score_r2g, e_score_g2r)
			logging.error('%s', e)
			sys.exit(1)
		if compound_score >= threshold_compound_score and reciprocal_score == threshold_reciprocal_score and e_score_r2g > threshold_e_score_r2g and e_score_g2r > threshold_e_score_g2r:
			outfile.write('\t'.join([gene_id, database_id_g2r]) + '\n')


def iter_columns_by_name(infile, col_names, delim=','):
	'''
	Will yield each line of infile, parsing only the columns whose headers were provided in col_names
	
	NOTE:
		- Assumes fist line contains column names
		- Will return lines in the order they appear in col_names
		- Will return an error if not all column headers were found.
		- Header names are case sensitive and must be complete word matches. 
	'''
	
	## Get first line from file. Assume it contains the column headers.
	header_line = infile.next().strip('\n')
	headers = header_line.strip('\n').split(delim)
	
	## Get the index of column names in infile.
	error_count = 0 # Keep track of errors
	headers_index = []
	for col_name in col_names:
		try:
			headers_index.append(headers.index(col_name))
		except ValueError:
			error_count += 1
			logging.error('Column name "%s" not found in header row of input file', col_name)
	## If we have encontered errors stop and print header row
	if error_count > 0:
		logging.error('Column names missing from infile. Stopping!')
		logging.error('Problem header row: %s', header_line)
		sys.exit(1)
	
	## Iterate over the file and return the columns of interest.
	for line in infile:
		line = line.strip('\n')
		if not line or line.startswith('#'):
			continue # Ignore blank or comment lines
		
		try:
			line_split = line.split(delim)
			yield [line_split[i] for i in headers_index]
		except IndexError:
			logging.error('Row found that is missing some columns!')
			logging.error('Problem row: %s', line)
			logging.error('Problem row split: %s', line_split)
			logging.error('Column indexes used: %s', headers_index)
			sys.exit(1)


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
