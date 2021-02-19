#vcfs = open("filtered_file_names.txt")
#fastas = open("fasta_file_names.txt")

###############################################
####makes a dictionary of chromosomes w/#######
####chromosome name as key and sequence########
####as value                           ########
###############################################

import sys


fasta_dict = {}
currentline = ""
fasta=open(sys.argv[1],'r')
#fasta = open('GCF_000599545.1_ASM59954v1_genomic.fna','r')

for line in fasta:
    if line.startswith('>'):
        line=line.rstrip('\n')
        if currentline != "" : fasta_dict[key] = currentline
        key = line[4:14]
#        print(key)
        currentline = ""
    else:
        line=line.rstrip('\n')
        currentline = currentline + line

fasta_dict[key] = currentline
currentline = ""
fasta.close()
bp = fasta_dict.get("CP003949.1")


###############################################
####CODON DICTIONARIES FOR TRANSLATION     ####
############################################### 

codon_dict = {'TTT': 'F', 'TTC':'F', 'TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y', 'TAC':'Y','TAA':'STOP','TAG':'STOP', 'TGT':'C','TGC':'C', 'TGA':'STOP','TGG':'W','CTT':'L','CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R', 'CGG':'R', 'ATT':'I', 'ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V', 'GTA':'V', 'GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A', 'GAT':'D', 'GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

start_codon_dict = {'TTG':'M*', 'CTG':'M*', 'ATT':'M*', 'ATC':'M*','ATA':'M*','ATG':'M*', 'GTG':'M*'}

###############################################
####makes the gff file into a list of lists####
############################################### 
gff = open(sys.argv[2],'r')
#gff = open('GCF_000599545.1_ASM59954v1_genomic.gff','r')
list_of_lists = []
for line in gff:
    if line[0] != "#":
        stripped_line = line.rstrip('\n')
        split_line = stripped_line.split('\t')
        if len(split_line) >= 8:
            list_of_lists.append(split_line)
gff.close

###############################################
####extracts snp data from vcf file        ####
###############################################
SNP_locations = open(sys.argv[3],'r')

#SNP_locations = open('locations.txt','r')
snps = []
for line in SNP_locations:
    stripped = line.strip('\n')
    snps.append(stripped)

BP_changes = open(sys.argv[4],'r')
#BP_changes = open('bp_changes.txt', 'r')
bp_change = []
for line in BP_changes:
    stripped = line.strip('\n')
    bp_change.append(stripped)

####Deprecated code that took hard inputs of <chromosome>_<bp.

#snps = ["CP003949.1_396714","CP003949.1_2370871","CP003949.1_3063788","CP003949.1_3306824","CP003949.1_4802458","CP003949.1_6359949","CP003949.1_8217888","CP003949.1_8326#930","CP003950.1_2231","CP003950.1_4656","CP003950.1_4872","CP003950.1_5254","CP003950.1_5488","CP003950.1_5699","CP003950.1_5818","CP003950.1_6107","CP003950.1_6182","CP#003950.1_6323","CP003950.1_6779","CP003950.1_6953","CP003950.1_7211","CP003950.1_13155","CP003950.1_13234","CP003950.1_13548","CP003950.1_13911","CP003950.1_19062","CP003#950.1_60604","CP003952.1_4420"]
    
#bp_change = ["C to T","T to G","C to A","G to A","G to A","G to T","C to T","C to A","C to T","G to A","T to G","A to G","G to C","G to C","T to C","G to A","G to A","G t#o A","T to C","T to C","G to A","A to G","C to T","T to C","G to C","A to G","A to G","T to C"]

z = -1

codon_output = "basepair, old codon, new codon, old AA, new AA, AA position \n"

#print(basepair,codon,new_codon_str,aa,new_aa,triple)

for i in snps:
    z+=1
    chromosome,basepair = i.split('_')
    for j in list_of_lists:
        if j[0] == "NZ_" + chromosome:
            if j[2] == 'gene' and int(basepair) >= int(j[3]) and int(basepair) <= int(j[4]):
                lower = int(j[3])
                upper = int(j[4])
                change = bp_change[z]
        #        z +=1
                old_bp,to,new_bp = change.split(' ')
#                print new_bp
                if j[6] == '+':
                    actual_lower = lower - 1
                    x = fasta_dict.get(chromosome)
                    gene = x[actual_lower:upper]
                    actual_basepair = int(basepair) - lower
                    triple = actual_basepair // 3
                    location = actual_basepair%3
                    codon = gene[triple*3:triple*3+3]
                    if actual_basepair <= 2:
                        aa = start_codon_dict.get(codon)
                    else:
                        aa = codon_dict.get(codon)
                    new_codon = list(codon)
                    new_codon[location] = new_bp
                    new_codon_str = ''.join(new_codon)
                    new_aa = codon_dict.get(new_codon_str)
                    triple += 1
                    print(basepair,codon,new_codon_str,aa,new_aa,triple)
                    
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
#                    print(gene_reversed)
                    triple = actual_basepair//3
                    location = actual_basepair%3
                    codon = gene_reversed[triple*3:triple*3+3]
                    if actual_basepair <= 2:
                        aa = start_codon_dict.get(codon)
                    else:
                        aa = codon_dict.get(codon)
                    
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
                    if actual_basepair <= 2:
                        new_aa = start_codon_dict.get(new_codon_str)
                    else:
                        new_aa = codon_dict.get(new_codon_str)
                    triple += 1
                    codon_output+= str(basepair)+" "+codon+" "+new_codon_str+" "+aa+" "+new_aa+" "+str(triple)+"\n"
                    #print(basepair,codon,new_codon_str,aa,new_aa,triple)
                
print(codon_output)
