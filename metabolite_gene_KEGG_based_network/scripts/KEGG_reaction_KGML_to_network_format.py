#!/usr/bin/env python2
DESCRIPTION = '''
Takes a KEGG reaction KGML file and reformats it for network visulization using Cytoscape.

NOTE:
	- KEGG reaction maps have compounds as nodes and reactions as edges.
	- Also adds 'gene' and 'ortholog' nodes.
	- Nodes without edges will have self edges so Cytoscape will plot them
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
~/GitHub/Construct_Networks/metabolite_gene_KEGG_based_network/scripts/KEGG_reaction_KGML_to_network_format.py			<subtype name="compound" value="XX"/>: The "compound" entry which links the two reactions.
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
	
	## Add edges for nodes without edges defined in KGML file.
	all_edges = []
	for x in edges:
		all_edges.extend(x)
	all_edges = set(all_edges)
	
	for n in nodes:
		n_id = n[0]
		if n_id not in all_edges:
			edges.append([n_id, n_id]) # Add edge to self
	
	## Write nodes to file
	nodes_file.write('node_id\tkegg_id\tname\ttype\tinfo\tlink\tx\ty\twidth\theight\tshape\n')
	for x in nodes:
		nodes_file.write('\t'.join([str(y) for y in x])+'\n')
	
	## Write edges to file
	edge_file.write('node_1\tnode_2\n')
	for x in edges:
		edge_file.write('\t'.join([str(y) for y in x])+'\n')




def parse_entry(entry):
	## Get node (compound, reaction or map) info
	# RETURN: [node_id, kegg_id, name, type, info, link, x, y, width, height, shape]
	
	####
	#### Get 'graphics' propoerties for each 'entry'
	####
	## Set defaults
	x_loc = 0
	y_loc = 0
	width = 46
	height = 17
	shape = 'circle'
	graphics_name = '-'
	name = '-'
	info = '-'
	## Access 'graphics' tag in 'entry'
	for x in entry:
		logging.debug('%s %s', x.tag, x.attrib) ## DEBUG
		if x.tag == 'graphics':
			## Check if we are dealing with explicit x and y positions or a list of coords.
			## If we are dealing with coords we should use the average as our x and y positions for plotting
			if 'x' in x.attrib and 'y' in x.attrib:
				x_loc = x.attrib['x']
				y_loc = x.attrib['y']
				width = x.attrib['width']
				height = x.attrib['height']
				if 'name' in x.attrib.keys():
					graphics_name = x.attrib['name']
			else:
				coords = x.attrib['coords']
				coords = coords.split(',')
				x_loc = sum([int(i) for i in coords[0::2]]) / len(coords[0::2])
				y_loc = sum([int(i) for i in coords[1::2]]) / len(coords[1::2])
				if 'name' in x.attrib.keys():
					graphics_name = x.attrib['name']
			shape = x.attrib['type']
	
	####
	#### Get info from kegg.jp link in 'entry'
	####
	## Stop if we dont have a link to work with. 
	if 'link' not in entry.attrib.keys():
		return [entry.attrib['id'], entry.attrib['name'], name, entry.attrib['type'], info, '-', x_loc, y_loc, width, height, shape]
	
	## Get info from link provided. 
	url = entry.attrib['link']
	html_text = requests.get(url).text
	soup = BeautifulSoup(html_text, 'html.parser')
	
	## Get table(s) from html
	tables = soup.findAll(lambda tag: tag.name=='table')
	## Check that we actually found some tables. If we didnt then we probabily have an incorrect URL or the entry doesn't have the info we are after.
	if not tables: # True if empty
		logging.info('No table found on KEGG website for entry "%s"', entry.attrib['name']) ## INFO
		logging.info('URL: "%s"', url) ## INFO
		return [entry.attrib['id'], entry.attrib['name'], '-', entry.attrib['type'], '-', '-', x_loc, y_loc, width, height, shape]
	table = tables[0]
	
	## Split master table into rows (whole webpage is a table so we have to access table within a cell of a table)
	master_rows = table.findAll(lambda tag: tag.name=='tr')
	rows = master_rows[0].findAll('td')[0].findAll('tr')[1].findAll('td')[0].findAll('tr')
	
	## type="compound": Get "Name" (first if multiple) and "Exact mass" from KEGG database.
	if entry.attrib['type'] == 'compound':
		# Iterate over rows and check if we have found the correct rows.
		for row in rows:
			# Get value from 1st and second columns (or return blank values if split failed)
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Name':
				name = row_value.split(';')[0]
			elif row_name == 'Exact mass':
				info = row_value
		## KEGG compound ID examples: cpd:C02226 (need to lstrip "cpd:")
		return [entry.attrib['id'], entry.attrib['name'].lstrip("cpd:"), name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="reaction": Get "Enzyme" from KEGG database.
	elif entry.attrib['type'] == 'reaction':
		name_list = []
		info_list = []
		for table in tables:
			master_rows = table.findAll(lambda tag: tag.name=='tr')
			try:
				rows = master_rows[0].findAll('td')[0].findAll('tr')[1].findAll('td')[0].findAll('tr')
			except IndexError:
				continue
			# Iterate over rows and check if we have found the correct rows (have to also strip \n and \xa0 characters from string).
			for row in rows:
				row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
				row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
				if row_name == 'Entry' and 'RClass' in row_value:
					break
				if row_name == 'Enzyme':
					name_list.extend(row_value.split()) # Remove large blocks of white space between multiple IDs
				if row_name == 'Orthology':
					tmp_value = "" if len(row.findAll('td')) == 0 else row.findAll('tr') # Get Orthology rows
					info_list.extend([x.get_text().strip().replace(u'\xa0', u' ').split(' ')[0] for x in tmp_value]) # Get row text, clean it, then split and take first work (which is 'K' ID)
		## KEGG reaction ID examples: rn:R05071 rc:RC00837 (split by space and need to lstrip "rn:")
		name = ';'.join(set(name_list))
		info = ';'.join(set(info_list))
		return [entry.attrib['id'], entry.attrib['name'].split(' ')[0].lstrip("rn:"), name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="map": Get "Name" (first if multiple) from KEGG database.
	elif entry.attrib['type'] == 'map':
		# Iterate over rows and check if we have found the correct rows.
		for row in rows:
			# Get value from 1st and second columns (or return blank values if split failed)
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Name':
				name = row_value.split(';')[0]
		## KEGG pathway ID example: path:rn00290 (need to lstrip "path:")
		return [entry.attrib['id'], entry.attrib['name'].lstrip('path:'), name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="gene": Get "Name" (first if multiple) from KEGG database.
	elif entry.attrib['type'] == 'gene':
		ids4info = []
		for table in tables:
			master_rows = table.findAll(lambda tag: tag.name=='tr')
			try:
				rows = master_rows[0].findAll('td')[0].findAll('tr')[1].findAll('td')[0].findAll('tr')
			except IndexError:
				continue
			for row in rows:
				row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
				if row_name == 'KO':
					row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
					ids4info.append(row_value.split(' ')[0])
		info = ';'.join(set(ids4info))
		return [entry.attrib['id'], entry.attrib['name'], graphics_name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## type="ortholog": Get "Name" (first if multiple) from KEGG database.
	elif entry.attrib['type'] == 'ortholog':
		# Iterate over rows and check if we have found the correct rows.
		for row in rows:
			# Get value from 1st and second columns (or return blank values if split failed)
			row_name = "" if len(row.findAll('th')) == 0 else row.findAll('th')[0].get_text().strip().replace(u'\xa0', u' ')
			if row_name == 'Definition':
				row_value = "" if len(row.findAll('td')) == 0 else row.findAll('td')[0].get_text().strip().replace(u'\xa0', u' ')
				name = row_value
		kegg_id = ';'.join([x.lstrip('ko:') for x in entry.attrib['name'].split(' ')])
		return [entry.attrib['id'], kegg_id, name, entry.attrib['type'], info, entry.attrib['link'], x_loc, y_loc, width, height, shape]
	
	## If we didnt account for something
	else:
		logging.info('Found entry "type" that we havent accounted for: %s', entry.attrib) ## INFO
	


def parse_relation(relation):
	## Get edge (relation) info
	## RETURNS: [node_1, node_2]
	##   node_1: "reaction" node id
	##   node_2: "compount" node id
	compound_node_id = 'undefined'
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
