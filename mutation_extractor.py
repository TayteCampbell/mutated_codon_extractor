##################################################
####makes a dictionary of chromosomes with########
####chromosome name as key and sequence   ########
####as value                              ########
##################################################

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description = """
    This script takes a reference nucleotide genome assembly, 
    a .gff file, and mutations from a .vcf file, and 
    returns the bp location, original codon, the mutated codon, the original
    amino acid, the mutated amino acid, and the numerical position of the amino acid
    within the protein as a pandas dataframe. Intergenic SNPs are not returned in the output.
    """
)

parser.add_argument(
    "-assembly", dest = "assembly", help =
    """This file should be the reference assembly that you used for variant calling""", 
    default=None, required=True
)

parser.add_argument(
    "-gff", dest = "gff", help =
    """This file should be a .gff file that lists chromosome and bp locations of annotated genes""",
    default=None, required=True
)

parser.add_argument(
    "-mut", dest = "mutations",help =
    """A comma-separated file with the basepair changes in the variant strain. For example:
    A,C
    T,G
    G,A
    """,
    default = None, required = True
)

parser.add_argument(
    "-loc", dest = "locations", help = 
    """This file should contain the chromosome bp locations of each SNP separated 
    by an underscore. For example:
    CP003950.1_1137
    CP003950.1_1448 
    """,
    default = None, required = True
)

inputs = parser.parse_args()

fasta_dict = {}
currentline = ""
with open(f'{inputs.assembly}', 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            line=line.rstrip('\n')
            if currentline != "" : fasta_dict[key] = currentline
            key = line[4:14]
            currentline = ""
        else:
            line=line.rstrip('\n')
            currentline = currentline + line

    fasta_dict[key] = currentline
    currentline = ""


######################################################################
#### Codon dictionaries for translation and translation function  ####
######################################################################

codon_dict = {'TTT': 'F', 'TTC':'F', 'TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y', 'TAC':'Y','TAA':'STOP','TAG':'STOP', 'TGT':'C','TGC':'C', 'TGA':'STOP','TGG':'W','CTT':'L','CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R', 'CGG':'R', 'ATT':'I', 'ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V', 'GTA':'V', 'GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A', 'GAT':'D', 'GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

start_codon_dict = {'TTG':'M*', 'CTG':'M*', 'ATT':'M*', 'ATC':'M*','ATA':'M*','ATG':'M*', 'GTG':'M*'}


def determine_aa_from_codon(bp_in_gene, codon_seq):
    if bp_in_gene <=2:
        amino_acid = start_codon_dict.get(codon_seq)
    else:
        amino_acid = codon_dict.get(codon_seq)
    return amino_acid


#################################################
#### makes the gff file into a list of lists ####
################################################# 
#gff = open(sys.argv[2],'r')
#gff = open('GCF_000599545.1_ASM59954v1_genomic.gff','r')

with open(f'{inputs.gff}','r') as gff:
    list_of_lists = []  
    for line in gff:
        if line[0] != "#":
            stripped_line = line.rstrip('\n')
            split_line = stripped_line.split('\t')
            if len(split_line) >= 8:
                list_of_lists.append(split_line)


##############################################################
#### extracts snp data from locations and mutations files ####
##############################################################


with open(f'{inputs.locations}','r') as SNP_locations:
    extracted_locations = []
    for line in SNP_locations:
        stripped = line.strip('\n')
        extracted_locations.append(stripped)
    

with open(f'{inputs.mutations}','r') as BP_changes:
    bp_change = []
    for line in BP_changes:
        stripped = line.strip('\n')
        bp_change.append(stripped)




codon_output = "basepair, old codon, new codon, old AA, new AA, AA position \n"

###Parse through the SNPs and find their locations in the genes using the .gff file

output_df = pd.DataFrame(columns = ['Chromosome', 'Basepair','old codon','new codon','old AA','new AA', 'AA position'])
z = -1
for i in extracted_locations:
    z += 1
    chromosome,basepair = i.split('_')
    ### Iterate through the .gff lines
    for j in list_of_lists:
        ### Find the correct chromosome
        if j[0] == "NZ_" + chromosome:
            ### Make sure we're looking at a gene line and if our bp locations is between the 
            ### start and stop locations then save the locations as lower and upper and get the 
            ### SNP bp change location
            if j[2] == 'gene' and int(basepair) >= int(j[3]) and int(basepair) <= int(j[4]):
                lower = int(j[3])
                upper = int(j[4])
                change = bp_change[z]
                old_bp,new_bp = change.split(',')
                ### Check column 6 in the .gff to see if on the forward or reverse strand
                if j[6] == '+':
                    actual_lower = lower - 1 ### change the lower value to align with python indexing values
                    x = fasta_dict.get(chromosome) ### get the DNA sequence of the chromosome
                    gene = x[actual_lower:upper] ### Slice the gene out of the chromosome
                    actual_basepair = int(basepair) - lower ### find the bp number for the gene
                    triple = actual_basepair // 3 ### get the codon number for the gene
                    location = actual_basepair%3 ### find where in the codon the mutation occurred
                    codon = gene[triple*3:triple*3+3] ### Slice the codon out of the gene
                    
                    ### get the amino acid from the codon
                    aa = determine_aa_from_codon(actual_basepair, codon)
                    
                    ### create the new codon with the mutation
                    new_codon = list(codon) 
                    new_codon[location] = new_bp
                    new_codon_str = ''.join(new_codon)

                    ### get the net amino acid
                    new_aa = determine_aa_from_codon(actual_basepair, new_codon_str)
                    triple += 1 #correct for python indexing
                    output_df.loc[len(output_df)] = [chromosome, basepair, codon, new_codon_str, aa, new_aa, triple]


                ### if the gene is on the reverse strand we need to get the reverse complement

                elif j[6] == '-':
                    actual_lower = lower-1
                    actual_basepair = upper-int(basepair)
                    x = fasta_dict.get(chromosome)
                    gene = x[actual_lower:upper]
                    reverse_gene = gene[::-1]
                    gene_reversed = ''
                    for i in reverse_gene:
                        if i == "A":
                            gene_reversed += "T"
                        if i == "T":
                            gene_reversed += "A"
                        if i == "C":
                            gene_reversed += "G"
                        if i == "G":
                            gene_reversed += "C"
                    triple = actual_basepair//3
                    location = actual_basepair%3
                    codon = gene_reversed[triple*3:triple*3+3]
                    
                    aa = determine_aa_from_codon(actual_basepair, codon)
                                        
                    new_bp_complement = ''
                    for i in new_bp:
                        if i == "A":
                            new_bp_complement += "T"
                        if i == "T":
                            new_bp_complement += "A"
                        if i == "C":
                            new_bp_complement += "G"
                        if i == "G":
                            new_bp_complement += "C"
                    new_codon = list(codon)
                    new_codon[location] = new_bp_complement
                    new_codon_str = ''.join(new_codon)
                    
                    new_aa = determine_aa_from_codon(actual_basepair, new_codon_str)
                    
                    triple += 1
                    output_df.loc[len(output_df)] = [chromosome, basepair, codon, new_codon_str, aa, new_aa, triple]


                

                
print(output_df)