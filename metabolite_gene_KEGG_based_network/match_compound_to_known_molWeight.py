#!/usr/bin/env python2
DESCRIPTION = '''
Takes a file with compound molecular weights and matches them to the closes molecular weight from a set of known mol weights.


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
	parser.add_argument('--unknown', metavar='unknonw.txt', 
		required=True, type=lambda x: File(x, 'r'), 
		help='Input [gzip] file with unknown mol weights (required)'
	)
	parser.add_argument('--known', metavar='known.txt',
		required=True, type=lambda x: File(x, 'r'),
		help='Input [gzip] file with known mol weights (required)'
	)
	parser.add_argument('-o', '--out', metavar='output.txt', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] file (default: stdout)'
	)
	parser.add_argument('--unknown_col', 
		required=True, metavar=1, type=int, 
		help='Column (1-based) to find unknown mol weights in (required)'
	)
	parser.add_argument('--known_col',
		required=True, metavar=1, type=int,
		help='Column (1-based) to find known mol weights in (required)'
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
	
	
	with args.unknown as unknown_file, args.known as known_file, args.out as out_file:
		## Dont forget to change from 0-based to 1-based index
		match_compound_to_known_molWeight(unknown_file, known_file, out_file, args.unknown_col-1, args.known_col-1)



def match_compound_to_known_molWeight(unknown_file, known_file, out_file, unknown_col, known_col, 
					max_diff=5.0, unknown_file_delim='\t', known_file_delim='\t'):
	
	known = {}
	for line in known_file:
		line = line.rstrip('\n')
		if not line or line.startswith('#'):
			continue
		
		line_split = line.split(known_file_delim)
		try:
			molWeight = float(line_split[known_col])
			known[molWeight] = line
		except IndexError:
			logging.info("ERROR: %s", line)
			logging.info("ERROR: --known_col %s out of range for --known", known_col+1)
			sys.exit(1)
		except ValueError:
			logging.info("ERROR: %s", line)
			logging.info("ERROR: Can't change '%s' to float in --known", line_split[known_col])
			sys.exit(1)
	known_molWeights = known.keys()
	
	for line in unknown_file:
		line = line.rstrip('\n')
		if not line or line.startswith('#'):
			continue
		
		line_split = line.split(unknown_file_delim)
		try:
			molWeight = float(line_split[unknown_col])
		except IndexError:
			logging.info("ERROR: %s", line)
			logging.info("ERROR: --unknown_col %s out of range for --unknown", unknown_col+1)
			sys.exit(1)
		except ValueError:
			logging.info("ERROR: %s", line)
			logging.info("ERROR: Can't change '%s' to float in --unknown", line_split[unknown_col])
			sys.exit(1)
		
		## Loop over known mol weights and find the closest
		smallest_diff = 999999999 # Set super high so any difference is smaller
		closest_known_MW = None
		
		logging.debug('Starting search with mol weight: %s', molWeight) ## DEBUG
		for mw in known_molWeights:
			diff = abs(molWeight-mw)
			if diff < smallest_diff and diff <= max_diff:
				smallest_diff = diff
				closest_known_MW = mw
				logging.debug('NEW: smallest diff (%s) and known mol weight (%s)', smallest_diff, closest_known_MW) ## DEBUG
		logging.debug('Smallest diff found is %s and belongs to known mol weight %s', smallest_diff, closest_known_MW) ## DEBUG
		if smallest_diff != 999999999:
			out_file.write(line+'\t'+known[closest_known_MW]+'\n')
		else:
			out_file.write(line+'\n')



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
