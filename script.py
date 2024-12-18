import argparse
from app.models import FileReader, Sequence, ORF, AminoAcidChain

parser = argparse.ArgumentParser(description="A simple script to find the SNPs between the longest ORF in the reference genome and a sub sequence in the position on a variance genome.")
parser.add_argument("ref", help="Name of .fna file containing a single reference genome.")
parser.add_argument("var", help="Name of .fna file containing a single variance genome.")
parser.add_argument("-v", "--verbose", help="Help increase program verbosity.", action="store_true")
args = parser.parse_args()

if args.verbose:
    pass

fr = FileReader()
reference = fr.extract_sequences(args.ref)
variance = fr.extract_sequences(args.var)

reference_genome = Sequence(reference[0][0], reference[0][1])
variance_genome = Sequence(variance[0][0], variance[0][1])

reference_orf = reference_genome.get_longest_orf()
variance_orf = variance_genome.get_longest_orf()
variance_aa = AminoAcidChain(orf=variance_orf)
#variance_orf = variance_genome.get_subsequence_as_ORF(reference_orf.index, reference_orf.length, 3)
print(reference_orf.find_snps(variance_orf))
print(reference_orf.description())
print(variance_orf.description())
print(variance_aa.description())