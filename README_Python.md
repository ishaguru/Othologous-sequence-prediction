Orthologous CRM prediction pipeline

This pipeline makes use of NCBI Blast and is written using Biopython. Required Python modules include pandas, matplotlib, pyplot, re, NCBIWWW, NCBIXML, and Entrez. The developed computational framework addresses the challenge of identifying orthologous sequences in closely related insect species

Functions:

run_blast(query_seq, database_name, entrez_query, nucl_penalty=-1, nucl_reward=1, word_size=7)
This function runs a BLAST search using the NCBIWWW.qblast() function from the Bio.Blast library. It takes the following parameters:
query_seq (str): path to a file containing the query sequence in FASTA format.
database_name (str): name of the database to search against.
entrez_query (str): an Entrez query to limit the search to a specific organism or group of organisms.
nucl_penalty (int): nucleotide mismatch penalty (default=-1).
nucl_reward (int): nucleotide match reward (default=1).
word_size (int): word size used for the search (default=7).
The function saves the BLAST search results in XML format to an output file.

remove_word_from_xml(xml_file, word, output_file) 
This function removes a given word from an XML file and saves the updated file. In this case, this function was added to remove the word 'CREATE_VIEW' as it did ot allow to parse the xml. It takes the following parameters:
xml_file (str): path to the XML file to modify.
word (str): the word to remove.
output_file (str): path to the output file.
The function removes all occurrences of the given word from the XML file and saves the updated file to the output path.

parse_results(xml_file, output_file, evalue_thresh)
This function parses the BLAST search results from an XML file and saves them to a tab-delimited file. It takes the following parameters:
xml_file (str): path to the XML file to parse.
output_file (str): path to the output file.
evalue_thresh (float): an E-value threshold for filtering the search results.
The function parses the XML file using the NCBIXML.parse() function and filters the results based on the given E-value threshold. The filtered results are saved to the output file in tab-delimited format with the following columns: ID, Start, Stop, Species, E-value, Score, and Align_length.

filter_score()
This function reads the output file from parse_results() and filters the results based on the score column. It also plots the score distribution using matplotlib.pyplot.scatter(). 

flank(bp, input_file, output_file)
This function retrieves sequences from the Entrez API based on the search results from parse_results(). It takes the following parameters:
bp (int): the number of base pairs to retrieve upstream and downstream of the sequence start and end positions.
input_file (str): path to the input file containing the search results in tab-delimited format.
output_file (str): path to the output file to save the retrieved sequences in FASTA format.
The function retrieves the sequences using the Entrez API efetch() function and saves them in FASTA format to the output file.
