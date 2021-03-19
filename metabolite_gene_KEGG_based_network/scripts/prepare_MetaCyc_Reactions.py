#!/usr/bin/env python2
DESCRIPTION = '''
Takes a big table with all reactions from Meatcyc and reformats ID mappings so they are nicer to work with. 

## Input file format (4 columns; expect header row).
MetaCyc [tab] KEGG [tab] RHEA [tab] EC-Number
...
...
FUCULOKIN-RXN [tab] <a href='http://www.genome.jp/dbget-bin/www_bget?rn:R03241'>R03241</a> [tab] <a href='http://www.rhea-db.org/reaction?id=12377'>12377</a> [tab] EC-2.7.1.51
RXN-18453  [tab]      [tab]     [tab]    EC-1.14.99.M4
RXN-11319  [tab] <a href='http://www.genome.jp/dbget-bin/www_bget?rn:R10246'>R10246</a> [tab] <a href='http://www.rhea-db.org/reaction?id=26364'>26364</a> [tab] EC-4.1.99.19
...
...

NOTE: Multiple IDs are seperated by ' // '

## Output:
MetaCyc_ID [tab] KEGG_Reaction_IDs [tab] RHEA_IDs [tab] EC_Numbers

KEGG_Reaction_IDs: Blank if no associated IDS, multiple IDs seperated by commas
RHEA_IDs: Blank if no associated IDS, multiple IDs seperated by commas
EC_Numbers: Blank if no associated IDS, multiple IDs seperated by commas
'''
import sys
import os
import argparse
import logging
import gzip
from bs4 import BeautifulSoup


## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-i', '--input', metavar='All_reactions_of_MetaCyc.txt.gz', 
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Input [gzip] MetaCyc ID mapping file (default: stdin)'
	)
	parser.add_argument('-o', '--out', metavar='MetaCyc_2_KEGG_Reaction_mapping.txt.gz', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] MetaCyc to KEGG Reaction mapping file (default: stdout)'
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
	
	
	with args.input as metaCycFile, args.out as outfile:
		process_MetaCyc_Reactions(metaCycFile, outfile)



def process_MetaCyc_Reactions(metaCycFile, outfile, writeDelim=','):
	'''
	Takes a file detailing all the reactions in MetaCyc and cleans up the ID mappings.
	
	Expected File format: MetaCyc [tab] KEGG [tab] RHEA [tab] EC-Number
	
	'''
	metaCycFile.next()
	for line in metaCycFile:
		line = line.strip('\n')
		if not line or line.startswith('#'):
			continue # Ignore blank or comment lines
		
		try:
			MetaCyc, KEGG_raw, RHEA_raw, EC_Number_raw = line.split('\t')
		except ValueError:
			logging.error("Filed to split line into 4 columns: %s", line)
			sys.exit(1)
		
		## Process IDS and write output
		##	- Seperate IDs split by " // "
		## 	- Also findAll each ID and iterate over each possible element.
		
		## KEGG IDs
		KEGG = []
		for IDs in KEGG_raw.split(' // '):
			for ID in BeautifulSoup(IDs, 'html.parser').findAll("a"):
				KEGG.append(ID.text)
		KEGG_string = writeDelim.join(KEGG)
		
		## RHEA IDs
		RHEA = []
		for IDs in RHEA_raw.split(' // '):
			for ID in BeautifulSoup(IDs, 'html.parser').findAll("a"):
				RHEA.append('RHEA:'+ID.text)
		RHEA_string = writeDelim.join(RHEA)
		
		## EC Numbers
		EC_Number = []
		for ID in EC_Number_raw.split(' // '):
			EC_Number.append(ID)
		EC_Number_string = writeDelim.join(EC_Number)
		
		## Write to output
		outfile.write('{}\t{}\t{}\t{}\n'.format(MetaCyc, KEGG_string, RHEA_string, EC_Number_string))



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
