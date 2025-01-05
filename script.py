import argparse

from snpdetector.models import ORF, SNP, AminoAcidChain, FileReader, Sequence

parser = argparse.ArgumentParser(
    description="A simple script to find the SNPs between the longest ORF in the reference genome and a sub sequence in the position on a variance genome."
)
parser.add_argument(
    "ref", help="Name of FASTA file containing a single reference genome."
)
parser.add_argument(
    "var", help="Name of FASTA file containing a single variance genome."
)
parser.add_argument(
    "-v",
    "--verbose",
    default=False,
    help="Help increase program verbosity.",
    action="store_true",
)
parser.add_argument(
    "-l",
    "--list",
    default=False,
    help="Include a list of all the SNPs in the output.",
    action="store_true",
)
parser.add_argument(
    "-o",
    "--output",
    default=False,
    help="Target filepath to write output as a text file.",
)
parser.add_argument(
    "-s",
    "--sameplace",
    default=False,
    help="Add this argument if you'd like to compare the longest ORF on the reference gene to a subsequence from the variance gene that is the same size, and starts from the same base on the same strand.",
    action="store_true",
)
args = parser.parse_args()

fr = FileReader()
reference = fr.extract_sequences(args.ref)
variance = fr.extract_sequences(args.var)

reference_genome = Sequence(reference[0][0], reference[0][1])
variance_genome = Sequence(variance[0][0], variance[0][1])

# Find ORFs
reference_orf = reference_genome.get_longest_orf()
if args.sameplace:
    variance_orf = variance_genome.get_subsequence_as_ORF(
        reference_orf.index, reference_orf.length, reference_orf.frame
    )
else:
    variance_orf = variance_genome.get_longest_orf()


# Create Amino Acid Chain objects
reference_aa = AminoAcidChain(orf=reference_orf)
variance_aa = AminoAcidChain(orf=variance_orf)

# Find and count the SNPs
snps = reference_orf.find_snps(variance_orf)
sub_count = 0
ins_count = 0
del_count = 0
hp_snps = []
ns_snps = []
s_snps = []
for snp in snps:
    if snp.severity == "High profile":
        hp_snps.append(snp)
    elif snp.severity == "Non-synonymous":
        ns_snps.append(snp)
    else:
        s_snps.append(snp)
    if snp.type == "Substitution":
        sub_count += 1
    elif snp.type == "Insertion":
        ins_count += 1
    else:
        del_count += 1

output = ""
output += "Reference Open Reading Frame:\n"
output += reference_orf.description() + "\n"
output += "Variance Open Reading Frame:\n"
output += variance_orf.description() + "\n"

# If the List argument has been specified, include a list of SNP descriptions to the output
if args.list:
    output += "High Profile Non-synonymous SNPs:\n"
    for snp in hp_snps:
        output += f"{snp.description(verbose=args.verbose)}\n"
    output += "Non-synonymous SNPs:\n"
    for snp in ns_snps:
        output += f"{snp.description(verbose=args.verbose)}\n"
    output += "Synonmous SNPs:\n"
    for snp in s_snps:
        output += f"{snp.description(verbose=args.verbose)}\n"

output += f"SNP severity count, Non-synonymous High profile - {len(hp_snps)}, Non-synonymous Low profile - {len(ns_snps)}, Synonymous - {len(s_snps)}\n"
output += f"SNP type count: Substitutions - {sub_count}, Insertions - {ins_count}, Deletions - {del_count}\n"

output += "Reference ORF as Amino Acid Chain:\n"
output += reference_aa.description() + "\n"
output += "Variance ORF as Amino Acid Chain:\n"
output += variance_aa.description() + "\n"

print(output)

if args.output:
    file = open(args.output, "w")
    file.write(output)
    file.close
