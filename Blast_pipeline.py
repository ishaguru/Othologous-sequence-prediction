
# Import required libraries
import pandas
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW, NCBIXML
import re
from Bio import Entrez

# Define function to run BLAST search and save results to a file
def run_blast(query_seq, database_name, entrez_query, nucl_penalty=-1, nucl_reward=1, word_size=7):
    """Performs a BLAST search and saves the results to a file."""
    with open(query_seq) as f:
        sequence_data = f.read()
        
# Perform BLAST search using NCBIWWW.qblast function and save results to a file

    result_handle = NCBIWWW.qblast(
        program="blastn",
        database=database_name,
        entrez_query=entrez_query,
        sequence=sequence_data,
        nucl_penalty=nucl_penalty,
        nucl_reward=nucl_reward,
        word_size=word_size
    )

    with open('results_2.xml', 'w') as f:
        f.write(result_handle.read())

# Define function to remove a word from an XML file and save the updated file    
def remove_word_from_xml(xml_file, word, output_file):
    """Removes a word from an XML file and saves the updated file."""
    with open(xml_file) as f:
        content = f.read()
    content = re.sub(word, '', content)
    with open(output_file, 'w') as f:
        f.write(content)

# Define function to parse BLAST results and save them to a file
def parse_results(xml_file, output_file, evalue_thresh):
    """Parses the BLAST results and saves them to a file."""
    with open(xml_file) as f:
        content = f.read()
    with open(output_file, 'w') as f:
        for record in NCBIXML.parse(open(xml_file)):
            if not record.alignments:
                continue
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < evalue_thresh:
                        hit_id = align.hit_id.split('|')
                        hit_def = align.hit_def.split(' ')
                        species = f"{hit_def[0]} {hit_def[1]}"
                        uid_new = hit_id[3]
                        f.write(f"{uid_new}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{species}\t{hsp.expect}\t{hsp.score}\t{hsp.align_length}\n")

# Define function to filter results based on score and plot the score distribution
def filter_score():
    # Read data from input file using pandas.read_csv function
    df = pandas.read_csv("output_11.txt", sep='\t',skiprows= [0],names=['ID', 'Start','Stop','Species','E-value','Score','Align_length','percent_identity'])
    data = df['Score']
    # Calculate mean, standard deviation, and upper limit based on standard deviation of the score distribution
    mean = df['Score'].mean()
    std = df['Score'].std()
    limit = mean + std
    # Determine minimum and maximum score values for plotting purposes
    min_value = min(data)
    max_value = max(data)
    # Plot the score distribution using matplotlib.pyplot.scatter and matplotlib
    plt.title("Score_values")
    plt.ylim(min_value - 100, max_value + 100)
    plt.scatter(x=df.index, y=df['Score'])
    plt.hlines(y= limit, xmin=0, xmax=len(data),colors='g')
    plt.show()
    
    print(limit)
    
    df_new = df[df['Score'] > limit]
    print(df_new)

    df_new.to_csv('final_op11.txt', sep="\t",index= False,header= False)

#This function reads the input file, extracts the ID, start and end positions of each sequence, fetches the corresponding sequences using the Entrez API, and writes the output to a file in FASTA format.
def flank(bp, input_file, output_file):
    uid = []
    start = []
    end = []
    with open(input_file) as f:
        rows = (line.split('\t') for line in f)
        for row in rows:
            uid.append(row[0])
            if int(row[1]) < int(row[2]):
                start.append(int(row[1]))
                end.append(int(row[2]))
            else:
                start.append(int(row[2]))
                end.append(int(row[1]))

    # Running efetch and returning the output in FASTA format
    Entrez.email = "ishaguru64@gmail.com"
    with open(output_file, "w") as bed:
        for i in range(len(uid)):
            try:
                handle_new = Entrez.efetch(db="nuccore", id=uid[i], rettype='fasta', retmode='txt', seq_start=start[i]-bp, seq_stop=end[i]+bp)
                print(handle_new.read(), file=bed)
            except Exception as e:
                print(f"Error fetching sequence {uid[i]}: {str(e)}")

    
#calling main function - example parameters  
if __name__ == "__main__":
    run_blast("agam_crm.txt" ,"refseq_genomic" ,"txid43816[ORGN]")
    remove_word_from_xml("results_2.xml", "CREATE_VIEW","file_2.xml")
    parse_results("file_2.xml", "output_2.txt", 0.05)
    filter_score()
    flank(10, "final_op2.txt", "agam_file")

