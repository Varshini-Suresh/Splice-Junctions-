#### DETERMINING LOCAITON OF INTRONS FROM SPLIT READS #### 

## IMPORT STATEMENTS ##
import sys
import logging
import re 

## SETTING UP THE LOGGER ##
logger = logging.getLogger()
logger.setLevel(logging.INFO) # For logger level INFO and higher messages 

ch = logging.StreamHandler() # setting up stream handler to avoid extra file outputs 
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

# testing the logger
logger.info("Testing that logging is working\n")


## WRITING FUNCTIONS ##
# Function to calculate location of junctions 
def find_junction(chr, pos, cigar):
    """
    Calculates start and end positions of junctions using regex
    Parameters:
      chr: chromosome number, where the read is aligned (str)
      pos: position on chromosome where alignment starts (str)
      cigar: string that describes how the read is aligned (str)
    Returns a nested dictionary with chr, start and end positions of the junctions 
    """
    junctions = {} 
    start = pos # setting start and end to position of the mapping 
    end = pos
    match = re.finditer(r'(\d+)([MINDS])', cigar)
    for m in match:
        num = int(m.group(1))
        char = m.group(2)

        # We only want to add matches (M) and deletions (D) to get starting position of the intron
        if char == 'M':
            start += num 
        if char == 'D':
            start += num
        # We need the skipped regions (N) to count end position of the intron 
        if char == 'N':
            end = start + num 
            key = str(chr) + '_' + str(start)+ '_' + str(end) # using chr, start and end for unique key for the outer dictionary
            junctions[key] = {'Chr': chr,'Start': start, 'End': end}
            start = end # setting end position of previous intron as start position for next possible N in the string
            end += num 
    return junctions

# Assert statements: for single and multiple Ns in the CIGAR string 
assert find_junction('TGME49_chrVIII', 300000, '63M146N37M') == {'TGME49_chrVIII_300063_300209': {'Chr': 'TGME49_chrVIII', 'Start': 300063, 'End': 300209}}
assert find_junction('TGME49_chrVIII', 654321, '45I348M75N1S48D59M578N10S') == {'TGME49_chrVIII_654669_654744': {'Chr': 'TGME49_chrVIII', 'Start': 654669, 'End': 654744}, 
                                                                                'TGME49_chrVIII_654851_655429': {'Chr': 'TGME49_chrVIII', 'Start': 654851, 'End': 655429}}

# Function to identify and count number of unique junctions 
def count_unique_junctions(junctions_dict_list):
    """
    Identifies and counts occurences of all unique junctions from a list of junction dictionaries
    Parameters: 
        junctions_dict_list: a list of nested dictionaries with details of all junctions from a given SAM file (list)
    Returns a nested dictionary containing details of only the unique junctions, along with number of each junction 
    present on the given list of junction dictionaries
    """
    unique_count = {}
    unique_junctions = {}
    count = 1

    # Creating two dictionaries - with unique junctions, and count 
    for dict in junctions_dict_list:
        for key, value in dict.items():
            if key not in unique_count:
                unique_count[key] = count
                unique_junctions[key] = value
            else: 
                unique_count[key] += count       
    
    # Adding count dictionary to the unique junctions dictionary
    for key, value in unique_junctions.items():
        data = unique_junctions[key]
        if unique_junctions.keys() == unique_count.keys():
            data['Count'] = unique_count[key]
    return unique_junctions
    
# Assert statements 
# Using junctions created from the previous function
a = find_junction('TGME49_chrVIII', 300000, '63M146N37M')
b = find_junction('TGME49_chrVIII', 654321, '45I348M75N1S48D59M578N10S') # Function should count for both junctions present in b

sample_junc_dict = [a, b, a, a, b]
c = count_unique_junctions(sample_junc_dict)

assert count_unique_junctions(sample_junc_dict) == {'TGME49_chrVIII_300063_300209': {'Chr': 'TGME49_chrVIII', 'Start': 300063, 'End': 300209, 'Count': 3}, 
                                                    'TGME49_chrVIII_654669_654744': {'Chr': 'TGME49_chrVIII', 'Start': 654669, 'End': 654744, 'Count': 2}, 
                                                    'TGME49_chrVIII_654851_655429': {'Chr': 'TGME49_chrVIII', 'Start': 654851, 'End': 655429, 'Count': 2}}

## TESTING INPUT FILES ## 
#logger.info("Please provide two files for input in the following order: the SAM file(.sam) and gene locations(.txt) \n")

# checking whether sys.argv has the right number of input files, if not issue logger error
if len(sys.argv) !=3: 
    sys.exit
    logger.error(f"Incorrect number of input files.\n" \
                 "Please provide two input files in the following order:\n" \
                 "the SAM file(.sam) and gene locations(.tsv) \n")
else: 
    logger.info(f"The files have been input correctly.\n")

# Assigning variables to the two input files 
sam_file = sys.argv[1]
gene_location_file = sys.argv[2]


## GETTING JUNCTIONS FROM THE GIVEN SAM FILE ## 
junction_list = [] # creating an empty list to store dict of junction locations

try: 
    with open(sam_file) as samfile: 

        # Extracting columns from the SAM file
        for row in samfile: 
            if row[0] != '@':
                row = row.rstrip().split('\t')
                rname = row[2]
                pos = int(row[3])
                cigar = row[5]
                nh_i_x = row[-1]

                # where read aligns only once i.e. where x = 1 in NH:i:x
                # len of string must be 6, to avoid other x ending 1, such as 11 or 21
                if len(nh_i_x)==6 and nh_i_x[-1] == '1': # here x is a string, not integer

                    # where the read is split i.e contains at least 1 N in CIGAR string
                    if 'N' in cigar:
                        junction = find_junction(rname, pos, cigar)
                        junction_list.append(junction)

except FileNotFoundError:
    sys.exit
    logger.error(f'File {sam_file} could not be opened. Please check and try again')

# With a list of junctions now created, the unique junctions can be filtered and counted
unique_junctions = count_unique_junctions(junction_list)


## GETTING GENE LOCATIONS FROM THE GIVEN TEXT FILE ## 
genes = {} #empty dict to store info on gene locations 

try: 
    with open(gene_location_file) as gene_location:
        header = next(gene_location)

        # Extracting columns from the tab file 
        for line in gene_location:
            line = line.rstrip().split('\t')
            geneID = line[0]
            transcriptID = line[1]
            location = re.split(r'[:..()]', line[2]) # Using regex to split string with multiple delimiters 
            chr = location[0]
            gene_start = int(location[1].replace(',', '')) # remove ',' from string
            gene_end = int(location[3].replace(',', ''))

            # Saving information in the genes dictionary 
            genes[geneID] = {"Chr": chr, "Start": gene_start, "End": gene_end}

except FileNotFoundError:
    sys.exit
    logger.error(f'File {gene_location_file} could not be opened. Please check and try again\n')


## WRITING LOCATIONS OF NUMBER OF UNIQUE JUNCTION PRESENT PER GENE, ONTO A TAB FILE ##

with open('Junctions.txt', 'w') as output: 

    # Assigning variables to gene location dictionary values
    for g_key, g_value in genes.items(): 
        gene_data = genes[g_key]
        gene_chr = gene_data['Chr']
        gene_start_pos = int(gene_data['Start'])
        gene_end_pos = int(gene_data['End'])

        # Assigning variables to junctions dictionary values 
        for j_key, J_value in unique_junctions.items():
            junc_data = unique_junctions[j_key]
            junc_chr = junc_data['Chr']
            junc_start = int(junc_data['Start'])
            junc_end = int(junc_data['End'])
            junc_count = junc_data['Count']

            # Checking if gene and junctions are on the same chromsome
            if junc_chr == gene_chr:
                # Junction location must be present within the boundaries of the gene 
                if junc_start >= gene_start_pos and junc_start <= gene_end_pos and junc_end <= gene_end_pos:
                    # If yes, write out the tab file 
                    output.write(f'{g_key}\t{junc_start}\t{junc_end}\t{junc_count}\n')
        output.write(f'\n') # Adding an empty row after the end of each gene 

# Telling user that the file ran with no errors 
logger.info('The file ran successfully. Output file is stored as Junctions.txt')

