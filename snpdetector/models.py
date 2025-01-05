nucleotide_complements = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

codon_table = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def is_nonsym_snp(three_mer_one, three_mer_two):
    """
    Determines if a Non-synonymous mutation has occured which results in different acids.
    """
    return True if codon_table[three_mer_one] != codon_table[three_mer_two] else False


def is_highprofile_snp(codon_one, codon_two):
    """
    A simple function that returns True if the difference between two codons is high profile.

    This function returns true if both codons are members of opposing acid groups.
    Otherwise returns False.
    """
    nonpolar_acids = ["G", "A", "V", "L", "I", "P", "F", "M", "W"]
    polar_acids = ["S", "T", "C", "N", "Q", "Y"]
    if (
        codon_table[codon_one] in nonpolar_acids
        and codon_table[codon_two] in polar_acids
    ):
        return True
    if (
        codon_table[codon_one] in polar_acids
        and codon_table[codon_two] in nonpolar_acids
    ):
        return True
    return False


class FileReader:
    """
    This class is a tool used to extract textual data from .fna files.

    The intended use of this class is to extract sequences from .fna files
    and return the sequences in a list of (header, sequence) tuples.
    """

    def extract_sequences(self, file_path):
        """
        Given a file path to a .fna file, trust that the file is stored in .fna format and extract all
        genomes from the sequence.

        The Filereader can extract more than one sequence from a .fna file.

        Returns a list of [[header][sequence],...]
        """
        file = open(file_path)
        lines = []
        while True:
            line = file.readline()
            if not line:
                break
            lines.append(line.strip())
        file.close()

        genomes = []
        sequence = ""
        while True:
            # Start from the end of the sequence until we find the header
            line = lines.pop(-1)
            if line.startswith(">"):
                genomes.append([line, sequence])
                sequence = ""
            else:
                sequence = line + sequence
            if not len(lines) > 0:
                break
        return genomes


class SNP:
    """
    A class to store data and functions for a Single Nucleotide Polymorphism(SNP)
    """

    def __init__(
        self, i_of_base, change_type, new_base="", severity="Synonamous", old_base=""
    ):
        self.index = i_of_base
        self.type = change_type
        self.new = new_base
        self.old = old_base
        self.severity = severity

    def description(self, verbose=False):
        """
        Returns a short or verbose description of the SNP object.
        """
        if verbose:
            description = f"SNP occured at base {self.index+1}(reference ORF). "
            if self.type == "Substitution":
                description += f"The base {self.old} in the reference was Substituted for {self.new} in the variance. "
            elif self.type == "Deletion":
                description += f"The base {self.old} in the reference was Deleted from the variance. "
            else:
                description += f"The base {self.new} was inserted in the variance. "
            if self.severity == "Synonymous":
                description += "The mutation was synonymous. "
            elif self.severity == "Non-synonymous":
                description += "The mutation was non-synonymous but low profile. "
            else:
                description += "The mutation was high profile. "
        else:
            description = (
                f"SNP location(reference): {self.index}, Type: {self.type}, Change: {self.old}->{self.new}, "
                f"Severity: {self.severity} "
            )
        return description


class ORF:
    """
    This is a model of an Open Reading Frame.
    An Open Reading Frame is a seuqnce of DNA that begins with a start codon
    and ends with a stop codon.
    """

    def __init__(self, sequence, start_nucleotide_index, orf_length, frame_translation):
        self.sequence = sequence
        self.index = start_nucleotide_index
        self.length = orf_length
        self.frame = frame_translation

    def description(self):
        """
        Returns a short or verbose description of the Open Reading Frame object.
        """
        description = (
            f"ORF starts at nucleotide {self.index+1} and contains {self.length} nucleotides.\n"
            f"ORF Sequence: {self.sequence}\nThe Sequence is on the "
        )
        description += "antisense strand.\n" if (self.frame) > 2 else "sense strand.\n"
        return description

    def is_longer(self, other_orf):
        return True if self.length >= other_orf.length else False

    def levenshtein_substitution_matrix(
        self, other_orf, insertion_penalty=1, deletion_penality=1, substituion_penalty=1
    ):
        """
        Create a substituion matrix using the Levenshtein algorithm and a reference sequence from the parent class and a
        variance sequence from an ORF provided as an argument.

        The penalties for an insertion, deletion, or substitution are modifiable depending on user need.

        Returns a substitution matrix.
        """
        reference = self.sequence
        variance = other_orf.sequence
        substitution_matrix = [
            [0 for _ in range(len(variance) + 1)] for _ in range(len(reference) + 1)
        ]

        for i in range(len(reference) + 1):
            for j in range(len(variance) + 1):
                if i == j == 0:
                    substitution_matrix[i][j] = 0
                elif 1 <= j <= len(variance) and i == 0:
                    substitution_matrix[i][j] = (
                        substitution_matrix[i][j - 1] + insertion_penalty
                    )
                elif 1 <= i <= len(reference) and j == 0:
                    substitution_matrix[i][j] = (
                        substitution_matrix[i - 1][j] + deletion_penality
                    )
                else:
                    if reference[i - 1] == variance[j - 1]:
                        substitution_matrix[i][j] = substitution_matrix[i - 1][j - 1]
                    else:
                        substitution_matrix[i][j] = min(
                            [
                                substitution_matrix[i][j - 1] + insertion_penalty,
                                substitution_matrix[i - 1][j] + deletion_penality,
                                substitution_matrix[i - 1][j - 1] + substituion_penalty,
                            ]
                        )
        return substitution_matrix

    def find_snps(self, other_orf):
        """
        This cumbersome function, bane of my existence, and foe of all who come after me, is surely beyond the comprehension of mortal ken.
        However I am requested to provide an explanation so I must do my best to explain what I have done.

        This function finds all the SNPs between the parent object and a compliment ORF object. It treats the parent class
        as the reference, and the compliment as the variant. To find the SNPs, it employs a helper function to produce a
        substitution matrix. It will greedily search for a path with the least SNPs on the substitution matrix. Can in theory, handle
        sequences of different length.

        When an SNP is found, it uses two helper functions to determine the severity of SNP before creating an SNP object
        and adding it to a list of SNPs.

        Returns a list of SNP objects.
        """
        sub_matrix = self.levenshtein_substitution_matrix(other_orf)
        snps = []
        i = j = 0
        # Whilst we have not reached the end of both sequences in our substitution matrix.
        while i < len(self.sequence) or j < len(other_orf.sequence):
            # The reference sequence still has bases left but the variance doesn't, the remaining bases in the reference were deleted.
            if j == len(other_orf.sequence):
                ref_offset = i % 3
                var_offset = j % 3
                ref_codon = self.sequence[i - ref_offset : i - ref_offset + 3]
                var_codon = other_orf.sequence[j - var_offset : j - var_offset + 3]
                snp = SNP(
                    i, "Deletion", old_base=self.sequence[i], severity="Non-synonymous"
                )
                snps.append(snp)
                i += 1
                continue
            # The reference sequence has bases left however the variance doesn't, the remaining bases in the variance were deleted.
            if i == len(self.sequence):
                ref_offset = i % 3
                var_offset = j % 3
                ref_codon = self.sequence[i - ref_offset : i - ref_offset + 3]
                var_codon = other_orf.sequence[j - var_offset : j - var_offset + 3]
                snp = SNP(
                    i,
                    "Insertion",
                    new_base=other_orf.sequence[j],
                    severity="Non-synonymous",
                )
                snps.append(snp)
                j += 1
                continue
            # Find the least "costly" move in our substituion matrix
            the_lowest_distance = min(
                sub_matrix[i + 1][j + 1], sub_matrix[i][j + 1], sub_matrix[i + 1][j]
            )
            # If the least costly move is a substitution, create a new SNP to record this substituion.
            if the_lowest_distance == sub_matrix[i + 1][j + 1]:
                if sub_matrix[i][j] != sub_matrix[i + 1][j + 1]:
                    ref_offset = i % 3
                    var_offset = j % 3
                    ref_codon = self.sequence[i - ref_offset : i - ref_offset + 3]
                    var_codon = other_orf.sequence[j - var_offset : j - var_offset + 3]
                    snp = SNP(
                        i,
                        "Substitution",
                        old_base=self.sequence[i],
                        new_base=other_orf.sequence[j],
                    )
                    if is_nonsym_snp(ref_codon, var_codon):
                        if is_highprofile_snp(ref_codon, var_codon):
                            snp.severity = "High profile"
                        else:
                            snp.severity = "Non-synonymous"
                    else:
                        snp.severity = "Synonymous"
                    snps.append(snp)
                i += 1
                j += 1
                continue
            # If the least costly move is an insertion, create a new SNP to record this insertion.
            if the_lowest_distance == sub_matrix[i][j + 1]:
                ref_offset = i % 3
                var_offset = j % 3
                ref_codon = self.sequence[i - ref_offset : i - ref_offset + 3]
                var_codon = other_orf.sequence[j - var_offset : j - var_offset + 3]
                snp = SNP(i, "Insertion", new_base=other_orf.sequence[j])
                if is_nonsym_snp(ref_codon, var_codon):
                    if is_highprofile_snp(ref_codon, var_codon):
                        snp.severity = "High profile"
                    else:
                        snp.severity = "Non-synonymous"
                else:
                    snp.severity = "Synonymous"
                snps.append(snp)
                j += 1
                continue
            # If the least costly move is a deletion, create a new SNP to record the deletion.
            if the_lowest_distance == sub_matrix[i + 1][j]:
                ref_offset = i % 3
                var_offset = j % 3
                ref_codon = self.sequence[i - ref_offset : i - ref_offset + 3]
                var_codon = other_orf.sequence[j - var_offset : j - var_offset + 3]
                snp = SNP(i, "Deletion", old_base=self.sequence[i])
                if is_nonsym_snp(ref_codon, var_codon):
                    if is_highprofile_snp(ref_codon, var_codon):
                        snp.severity = "High profile"
                    else:
                        snp.severity = "Non-synonymous"
                else:
                    snp.severity = "Synonymous"
                snps.append(snp)
                i += 1
                continue
        return snps

    def get_amino_acid(self):
        """
        A small function used to Create an amino acid sequence from the sequence in this ORF
        """
        aa_sequence = ""
        for i in range(0, len(self.sequence), 3):
            codon = codon_table[self.sequence[i : i + 3]]
            aa_sequence += codon
        return aa_sequence


class AminoAcidChain:
    """
    A class for storing amino acid chain sequences.
    """

    def __init__(self, acid_sequence="", orf=None):
        """
        This function will initialise an AminoAcidObject with the provided amino acid chain or ORF if provided with an orf.
        """
        if orf:
            self.sequence = orf.get_amino_acid()
        else:
            self.sequence = acid_sequence

    def description(self):
        """
        This class will generate a short or verbose description of the data contained in this AminoAcidChain object.
        """
        description = (
            f"This amino acid chain is {len(self.sequence)} acids long.\n"
            f"Sequence: {self.sequence}"
        )
        return description


class Sequence:
    """
    This is a model of DNA and mRNA sequences. This model wraps up functions
    useful for processing these sequences.
    """

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.anti_sense = ""
        for char in reversed(self.sequence):
            self.anti_sense += nucleotide_complements[char]

    def get_sequence_length(self):
        """
        Returns the length of the sequence.
        """
        return len(self.sequence)

    def get_longest_orf(self):
        """
        This function uses six frame translation to search for the longest Open Reading Frame.
        The a reverse compliment or anti-sense strand is produced from the sequence stored in the object
        and both are searched, three times each, to account for any frame shift that needs to be done
        to find the longest ORF.

        Returns an ORF object containing the Open Reading Frame, it's index, length, and frame shift used to find the ORF.
        """
        longest_orf = ORF("", 0, 0, 0)
        for translation in range(6):
            if translation > 2:
                search_area = self.anti_sense
            else:
                search_area = self.sequence
            for i in range(translation % 3, len(search_area), 3):
                head_codon = search_area[i : i + 3]
                if len(head_codon) > 2 and codon_table[head_codon] == "M":
                    for j in range(i + 3, len(search_area) - 2, 3):
                        tail_codon = search_area[j : j + 3]
                        if codon_table[tail_codon] == "*":
                            orf = ORF(search_area[i : j + 3], i, j - i + 3, translation)
                            if orf.is_longer(longest_orf):
                                longest_orf = orf
                            break
        return longest_orf

    def get_subsequence_as_ORF(self, start_index, length, translation):
        """
        This function is used to return the sequence or a subsequence as an ORF object.
        The index(starting at 0) of the base, the length of the subsequence, and the translation
        are used to extract the subsequence.

        A translation/frame shift 0-2 means the sub sequence will be taken from the sequence/sense strand.
        A translation/frameshift of 3-5 means the sub sequence will be taken from a reverse compliment/antisense strand
        of the stored sequence.

        Return an ORF object containing the subsequence.
        """
        if translation > 2:
            search_area = self.anti_sense
        else:
            search_area = self.sequence
        return ORF(
            search_area[start_index : start_index + length],
            start_index,
            length,
            translation,
        )
