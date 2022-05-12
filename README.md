# Mutated Codon Extractor
This is a personal script I wrote that takes point mutations extracted from a .vcf file and identifies the codons altered by each mutation.


Disclaimer: this script was written for a specific reference genome and .gff file for Rhodococcus opacus. Minor modifications may be needed to optimize this code for general purpose.

To use: 
1. Extract the mutation data from your .vcf file and put the chromosome and bp location in a locations file using the format "chromosome_bp". 
2. Extract the bp change from your .vcf file and put the old bp and new bp in the bp_changes file using the format "old bp,new bp".
3. Run the program with the .fna reference file, .gff reference file, locations file, and bp_changes file. For example:

mutation_extractor.py -assembly GCF_000599545.1_ASM59954v1_genomic.fna -gff GCF_000599545.1_ASM59954v1_genomic.gff -loc locations.txt -mut bp_changes.txt

4. The output will give you the chromosome, basepair location number, the old codon, the new codon (with the mutation), the old amino acid, the new amino acid, and the amino acid position in its respective protein.
 
Note: intergenic mutations will not show up in the output since they do not make any changes to genes.
  
