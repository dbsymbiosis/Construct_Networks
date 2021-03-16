#!/usr/bin/env python2
DESCRIPTION = '''
Takes a KEGG reaction KGML file and reformats it for network visulization using Cytoscape.

NOTE:
	- KEGG reaction maps have compounds as nodes and reactions as edges.
'''
import sys
import os
import argparse
import logging
import gzip
import requests
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET

## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-x', '--kgml', metavar='R01215.kgml',
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='KEGG reaction KGML file (default: stdin)'
	)
	parser.add_argument('-e', '--edges', metavar='network.edges.txt',
		required=True, default=sys.stdout, type=lambda x: File(x, 'w'),
		help='Output [gzip] network edges file (required)'
	)
	parser.add_argument('-n', '--nodes', metavar='network.nodes.txt',
		required=True, default=sys.stdout, type=lambda x: File(x, 'w'),
		help='Output [gzip] network node info file (required)'
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
	
	
	with args.kgml as kgml_file, args.edges as edge_file, args.nodes as nodes_file:
		KEGG_KGML_2_Cytoscape_network(kgml_file, edge_file, nodes_file)



def KEGG_KGML_2_Cytoscape_network(kgml_file, edge_file, nodes_file):
	'''
	entry: Info about a node or edge in the network
		type="reaction": KEGG reaction info (e.g. a enzymatic reaction)
		type="compound": KEGG compound info
		type="map":      Link to another KEGG map
	
	relation: Info about edges in the network (joining two reactions together)
		type="ECrel": Edge joining two "reaction" entries together
			<subtype name="compound" value="XX"/>: The "compound" entry which links the two reactions.
		type="maplink": Link to another map
	
	reaction: Info about the direction/substrate(s)/product(s) of a reaction
		type="irreversible"
		type="reversible"
	
	
	nodes:
		"compound" OR "reaction" OR "map" type entries
	edges:
		"relation" and "reaction" type entries (Need both as some compounds arent covered by "relation" tags only "reaction" tags)
	'''
	tree = ET.parse(kgml_file)
	root = tree.getroot()
	title = root.attrib['title']
	
	nodes = [] # all 'entry' tags
	edges = [] # build from 'relation' and 'reaction' tags
	
	for child in root:
		logging.debug('%s %s', child.tag, child.attrib) ## DEBUG
		if child.tag == 'entry':
			nodes.append(parse_entry(child))
			logging.debug("Entry parsed: %s", nodes[-1])
		elif child.tag == 'relation':
			for edge in parse_relation(child):
				## Add edge from 'relation' if not in list
				## [multiple 'relation' edges can pass through the same compound so we can end up with duplicate edges]
				if edge not in edges:
					edges.append(edge)
					logging.debug("Relation parsed and new edge added: %s", edge)
		elif child.tag == 'reaction':
			for edge in parse_reaction(child):
				## Add edge from 'reaction' if not in list (i.e. compound enters or exits the network)
				## [will not be in 'relation' tags becuase it is not between two reactions]
				if edge not in edges:
					edges.append(edge)
					logging.debug("Reaction parsed and new edge added: %s", edge)
		else:
			logging.info('Found tag in xml file that we havent accounted for: %s', child.attrib) ## INFO
		
	## Write nodes to file
	nodes_file.write('node_id\tkegg_id\tname\ttype\tinfo\tlink\tx\ty\twidth\theight\tshape\n')
	for x in nodes:
		nodes_file.write('\t'.join(x)+'\n')
	
	## Write edges to file
	edge_file.write('node_1\tnode_2\n')
	for x in edges:
		edge_file.write('\t'.join(x)+'\n')




def parse_entry(entry):
	## Get node (compound, reaction or map) info
	# RETURN: [node_id, kegg_id, name, type, info, link, x, y, width, height, shape]
	
	####
	#### Get 'graphics' propoerties for each 'entry'
	####
	## Set defaults
	x_loc = 0
	y_loc = 0
	width = 0
	height = 0
	shape = 'circle'
	## Access 'graphics' tag in 'entry'
	for x in entry:
		if x.tag == 'graphics':
			x_loc = x.attrib['x']
			y_loc = x.attrib['y']
			width = x.attrib['width']
			height = x.attrib['height']
			shape = x.attrib['type']
	
	####
	#### Get info from jegg.jp link in 'entry'
	####
	url = entry.attrib['link']
	html_text = requests.get(url).text
	soup = BeautifulSoup(html_text, 'html.parser')
	
	## Get table(s) from html
	table = soup.find(lambda tag: tag.name=='table')
	
	## Check that we actually found a table. If we didnt then we probabily used an incorrect URL
	if table is None:
		logging.error('No table found. URL found in XML file might be incorrect?') ## ERROR
		logging.error('URL found: %s', url) ## ERROR
		sys.exit(1)
	
	## Split master table into rows (whole webpage is a table so we have to access table within a cell of a table)
	master_rows = table.findAll(lambda tag: tag.name=='tr')
	rows = master_rows[0].findAll('td')[0].findAll('tr')[1].findAll('td')[0].findAll('tr')
	
	## type="compound": Get "Name" (first if multiple) and "Exact mass" from KEGG database.
	if entry.attrib['type'] == 'compound':
		info = 'NA'
		name = 'NA'
		# Iterate over rows and check if we have found the correct rows.
		for row in rows:
			# Get value from 1st and second columns (or return blank values if split failed)
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Name':
				name = row_value.split(';')[0]
			elif row_name == 'Exact mass':
				info = row_value
		return [entry.attrib['id'], entry.attrib['name'], name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="reaction": Get "Enzyme" from KEGG database.
	elif entry.attrib['type'] == 'reaction':
		info = 'NA'
		name = 'NA'
		# Iterate over rows and check if we have found the correct rows (have to also strip \n and \xa0 characters from string).
		for row in rows:
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Enzyme':
				name = ';'.join(row_value.split()) # Remove large blocks of white space between multiple IDs
		return [entry.attrib['id'], entry.attrib['name'], name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="map": Get "Name" (first if multiple) from KEGG database.
	elif entry.attrib['type'] == 'map':
		info = "NA"
		name = 'NA'
		# Iterate over rows and check if we have found the correct rows.
		for row in rows:
			# Get value from 1st and second columns (or return blank values if split failed)
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Name':
				name = row_value.split(';')[0]

		return [entry.attrib['id'], entry.attrib['name'], name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## If we didnt account for something
	else:
		logging.info('Found entry "type" that we havent accounted for: %s', entry.attrib) ## INFO
	


def parse_relation(relation):
	## Get edge (relation) info
	## RETURNS: [node_1, node_2]
	##   node_1: "reaction" node id
	##   node_2: "compount" node id
	for x in relation:
		if x.tag == 'subtype':
			compound_node_id = x.attrib['value']
	return [[relation.attrib['entry1'], compound_node_id],
		[relation.attrib['entry2'], compound_node_id]]


def parse_reaction(reaction):
	## Get edge info from reaction tags
	## 	'relation' tags give reaction<->compound<->reaction links
	##	'reaction' tags can give reaction<->compound links where the 
	##	 compound enters of leaves the network without a second reaction for 
	##	 the 'relation' tag to describe
	## RETURNS: [node_1, node_2]
	ret = []
	for x in reaction:
		if x.tag == 'substrate' or x.tag == 'product':
			ret.append([reaction.attrib['id'], x.attrib['id']])
	return ret



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
