# SNP Detector

## Use

    To use the script, place the script and the package in the same snpdetector file in the same directory. The snpdetector can only read sequences in the FASTA file format. The sequences use the standard dna codon table, therefore to use the script to examine mRNA sequences, the sequence has to be preprocessed to replace Uracil(U) bases with Thymine(T) bases.

    The files containing the reference and variant sequence are given as positional arguments. Ensure you are running the script from the directory with the files or use the absolute filepath destination to your two files.

    python ./script reference.fa variance.fa

    ## Optional arguments

    To see the complete list of optional arguments, run python ./script -h. 

    -v specifies the verbosity of descriptions and should be included with -l to produce a list of verbose descriptions of SNPs discovered

    -l can be added to include a list of descriptions for all the SNPs discovered, including the location of the SNP, it's type, the change, and so forth.

    -o can be added to specify a filepath to save the output to.

    -s This argument can be added if you'd like to assume the longest ORF on the variant gene starts from the same place on the same strand and is of the same length as the longest ORF on the reference genome.
