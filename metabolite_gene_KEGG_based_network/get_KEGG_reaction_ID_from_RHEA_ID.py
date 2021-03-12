#!/usr/bin/env python2
DESCRIPTION = '''
Takes a RHEA ID and scrapes the rhea-db.org website for the associated KEGG reaction ID

NOTE:
	- Assume all RHEA IDs have a **SINGLE** KEGG ID in a single cell.
	- Does **NOT** consider the direction of the reaction. 
		RHEA ID might be for a reaction the goes in a dirrection that is different to the KEGG reaction returned.
	- **CAN NOT** fetch KEGG reaction ID for non RHEA IDs (e.g. 5-NUCLEOTID-RXN). Will print en error and move on.
'''
import sys
import os
import argparse
import logging
import gzip
import requests
from bs4 import BeautifulSoup

## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-r', '--rhea', metavar='rhea_ids.txt',
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='RHEA ID to get KEGG reaction ID for (default: stdin)'
	)
	parser.add_argument('-o', '--out', metavar='output.txt',
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'),
	help='Output [gzip] file (default: stdout)'
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
	
	
	with args.rhea as infile, args.out as outfile:
		scrape_rheaDB(infile, outfile)



def scrape_rheaDB(RHEA_IDs, outfile):
	stop_on_invalid_url = False
	
	for RHEA_ID in RHEA_IDs:
		RHEA_ID = RHEA_ID.strip()
		if not RHEA_ID or RHEA_ID.startswith('#'):
			continue
		
		if ':' in RHEA_ID:
			rhea_url = "https://www.rhea-db.org/rhea/" + RHEA_ID.split(':')[1]
		else:
			rhea_url = "https://www.rhea-db.org/rhea/" + RHEA_ID
		logging.info('Scraping reaction info for %s from %s', RHEA_ID, rhea_url) ## DEBUG
		
		## Get "Links to other resources" table
		html_text = requests.get(rhea_url).text
		soup = BeautifulSoup(html_text, 'html.parser')
		table = soup.find(lambda tag: tag.name=='table' and tag.has_attr('id') and tag['id']=="otherresources") 
		
		## Check that we actually found a table. If we didnt then we probabily used an incorrect RHEA ID
		if table is None:
			logging.error('No "Links to other resources" table found. Maybe "%s" is not a valid RHEA ID?', RHEA_ID) ## ERROR
			if stop_on_invalid_url:
				sys.exit(1)
			else:
				continue
			
		
		## Split table into rows
		rows = table.findAll(lambda tag: tag.name=='tr')
		
		## Get first/header row and extract RHEA ids from each column
		col_names = []
		header_cells = rows[0].findAll(lambda tag: tag.name=='th')
		for cell in header_cells[1:]:
			# Get either <b> (for other related RHEA IDs) or <span> (for current RHEA ID 
			# [not a hyperlink which is why it uses <span> and not <b>])
			col_names.append(cell.findAll(lambda tag: tag.name=='b' or tag.name=='span')[0].get_text())
		logging.debug('"Links to other resources" table column names: %s', col_names) ## DEBUG
		
		## Iterate over rows to find 'KEGG' row
		for row in rows[1:]:
			cells = row.findAll(lambda tag: tag.name=='td')
			row_name = cells[0].findAll(lambda tag: tag.name=='b')[0].get_text()
			logging.debug('Row name: %s', row_name) ## DEBUG
			if row_name == "KEGG":
				## Find non-empty cells and write KEGG reaction ID to output
				for cell in cells[1:]:
					c = cell.get_text().strip()
					if c:
						outfile.write(RHEA_ID+"\t"+c+"\n")


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
