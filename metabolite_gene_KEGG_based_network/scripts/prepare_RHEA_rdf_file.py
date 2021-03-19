#!/usr/bin/env python2
DESCRIPTION = '''
Takes a big table with all reactions from RHEA and reformats ID mappings so they are nicer to work with. 

# Output:
RHEA_ID [tab] KEGG_Reaction_IDs [tab] MeatCyc_IDs [tab] EC_Numbers

KEGG_Reaction_IDs: Blank if no associated IDS, multiple IDs seperated by commas
MetaCyc_IDs: Blank if no associated IDS, multiple IDs seperated by commas
EC_Numbers: Blank if no associated IDS, multiple IDs seperated by commas
'''
import sys
import os
import argparse
import logging
import gzip
from lxml import etree


## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-r', '--rdf', metavar='rhea.rdf.gz', 
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Input [gzip] RHEA ID mapping file (default: stdin)'
	)
	parser.add_argument('-o', '--out', metavar='RHEA_2_KEGG_Reaction_mapping.txt.gz', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] RHEA to KEGG Reaction mapping file (default: stdout)'
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
	
	
	with args.rdf as rdffile, args.out as outfile:
		process_RHEA_rdf_file(rdffile, outfile)



def process_RHEA_rdf_file(rdffile, outfile):
	'''
	rdffile: RHEA RDF file from https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz
		Assumes undefined RHEA ID occures before directional IDs.
	
	outfile: Mapping between RHEA IDs, MetaCyc IDs and KEGG Reaction IDs
	'''
	RHEA_ID_annots = {} # {RHEA_ID_1:<annot class 001>, RHEA_ID_2:<annot class 001>, RHEA_ID_3:<annot class 002>}
	RHEA_DB_URL = 'http://rdf.rhea-db.org/'
	IDENTIFIERS_URL = 'http://identifiers.org/'
	IDENTIFIERS_DELIM = '/'
	EC_URL = 'http://purl.uniprot.org/'
	
	tree = etree.parse(rdffile)
	root = tree.getroot()
	nsmap = root.nsmap
	
	for child in root.findall('rdf:Description', nsmap):
		logging.debug('Processing element:\n%s', etree.tostring(child, pretty_print=True)) ## DEBUG
		
		## Get name (RHEA ID) of element.
		element_name = fetch_element_attrib(child, '{%s}about' % nsmap['rdf'], to_lstrip='http://rdf.rhea-db.org/')
		
		## Get the type of (subClassOf) reaction we are dealing with.
		## "Reaction":
		##     	- top level, will have elements describing other directional reactions.
		##     	- check for annotations
		## "BidirectionalReaction" and "DirectionalReaction":
		##     	- check for annotations
		subClassOf = fetch_attrib_from_elements(child, '{%s}subClassOf' % nsmap['rdfs'], '{%s}resource' % nsmap['rdf'], RHEA_DB_URL)
		if "Reaction" in subClassOf:
			## Associate this element's ID and sub-elements ID with same instance of annotation class
			RHEA_ID_annots[element_name] = Reaction_group_annotations()
			for react in fetch_attrib_from_elements(child, '{%s}directionalReaction' % nsmap['rh'], '{%s}resource' % nsmap['rdf'], RHEA_DB_URL):
				RHEA_ID_annots[react] = RHEA_ID_annots[element_name]
			for react in fetch_attrib_from_elements(child, '{%s}bidirectionalReaction' % nsmap['rh'], '{%s}resource' % nsmap['rdf'], RHEA_DB_URL):
				RHEA_ID_annots[react] = RHEA_ID_annots[element_name]
			
		elif "BidirectionalReaction" in subClassOf or "DirectionalReaction" in subClassOf:
			## Do nothing special, IDs already should be in "RHEA_ID_annots"; move onto next part where we find KEGG and MetaCyc annotations.
			pass
		else:
			## Found another type of element that we dont need (e.g. compound information)
			## Dont continue to the next part where we try and find links to external databases
			continue
			
		## Get links to other resources (if associated with this element)
		## Loop over each link and check if it a KEGG link, a MetaCyc link, or a EC link.
		attrib_list = []
		attrib_list.extend(fetch_attrib_from_elements(child, '{%s}seeAlso' % nsmap['rdfs'], '{%s}resource' % nsmap['rdf'], IDENTIFIERS_URL))
		attrib_list.extend(fetch_attrib_from_elements(child, '{%s}ec' % nsmap['rh'], '{%s}resource' % nsmap['rdf'], EC_URL))
		for attrib in attrib_list:
			## Expected attrib:
			##      "kegg.reaction/R02511"
			##     	"biocyc/METACYC:CYSTEAMINE-DIOXYGENASE-RXN"
			##      "biocyc/ECOCYC:1.5.8.2-RXN"
			try:
				database_name, database_id = attrib.split(IDENTIFIERS_DELIM)
				## Update annot instance with new annotation IDs
				if database_name == 'kegg.reaction':
					RHEA_ID_annots[element_name].add_KEGG(database_id)
				elif database_name == 'biocyc' and 'METACYC:' in database_id:
					## We only care about the 'METACYC' IDs for now.
					RHEA_ID_annots[element_name].add_MetaCyc(database_id.lstrip('METACYC:'))
				elif database_name == 'enzyme':
					RHEA_ID_annots[element_name].add_EC(database_id)
				else:
					pass
			except ValueError:
				logging.error('Could not split link "%s" to other resources corrrectly\n%s',
					attrib, etree.tostring(child, pretty_print=True)) ## ERROR
				sys.exit(1)
			except KeyError:
				logging.error('ID "%s" was not seen previously. This should not have happened.\n%s',
					attrib, etree.tostring(child, pretty_print=True)) ## ERROR
				sys.exit(1)
		
	## Write annotations to file.
	for annot_key, annot_value in RHEA_ID_annots.iteritems():
		outfile.write("RHEA:" + annot_key + '\t' + annot_value.writeKEGG() + '\t' + annot_value.writeMetaCyc() + '\t' + annot_value.writeEC() + '\n')



def fetch_attrib_from_elements(element, element_name, attrib_key, to_lstrip='', to_rstrip='', errorEmpty=False):
	'''
	Takes an element and finds all sub elements using element_name:
		- Throw error if nothing found and errorEmpty=True.
	For each sub-element try and fetch value of attrib_key:
		- Throw error if sub-element does not have attrib_key and errorEmpty=True.
	'''
	attrib_list = []
	## findall: will return a empty list of nothing found; return an error if empty list.
	sub_elements = element.findall(element_name)
	if len(sub_elements) == 0 and errorEmpty:
		logging.error('Could not find sub-elements with name "%s" in\n%s',
			element_name, etree.tostring(element, pretty_print=True)) ## ERROR
		sys.exit(1)
	for child in sub_elements:
		attrib = fetch_element_attrib(child, attrib_key, to_lstrip)
		attrib_list.append(attrib)
		logging.debug('Attrib found: %s', attrib) ## DEBUG
	return attrib_list



def fetch_element_attrib(element, attrib_key, to_lstrip='', to_rstrip=''):
	'''
	Fetch attrib value using attrib_key from element.
	Catch KeyError (key missing from attrib) and strip leading and trailing characters.
	'''
	try:
		attrib = element.attrib[attrib_key]
		attrib = attrib.lstrip(to_lstrip).rstrip(to_rstrip)
	except KeyError:
		logging.error('Cound not find "%s" attribute in element:\n%s',
			attrib_key, etree.tostring(element, pretty_print=True)) ## ERROR
	return attrib



class Reaction_group_annotations():
	'''
	Class to store the annotations from a group of RHEA IDS.
	
	- Specifically, 4 RHEA IDs represent the same reaction because they represent the reaction
		in a different direction/state (bidirectional, left-to-right, right-to-left, and undefined).
	- Main point of this class is that we can create a dict with RHEA IDs and associate all IDs with the
		same instance of this class. That way when we easily access the shared annotation for a set of
		IDs froma dict.
	'''
	def __init__(self):
		self.KEGG = [] # Store KEGG Reaction IDs
		self.MetaCyc = [] # Store Metacyc IDs
		self.EC = [] # Store EC numbers
	def add_KEGG(self, annot):
		if annot not in self.KEGG:
			self.KEGG.append(annot)
	def add_MetaCyc(self, annot):
		if annot not in self.MetaCyc:
			self.MetaCyc.append(annot)
	def add_EC(self, annot):
		if annot not in self.EC:
			self.EC.append(annot)
	def writeKEGG(self, delim=','):
		return delim.join(self.KEGG)
	def writeMetaCyc(self, delim=','):
		return delim.join(self.MetaCyc)
	def writeEC(self, delim=','):
		return delim.join(self.EC)



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
