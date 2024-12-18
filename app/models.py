nucleotide_complements = {
    'A': 'T', 'T': 'A',
    
    'G': 'C', 'C': 'G',
}

codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

class FileReader:
    """
    This class is a tool used to extract textual data from .fa files.
    
    The intended use of this class is to extract sequences from .fa files
    and return the sequences in a list of (header, sequence) tuples.
    """
    def extract_sequences(self, file_path):
        """_summary_

        Args:
            file_path (String): Path to .fa file containing sequence data

        Returns:
            (String)
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
            if line.startswith('>'):
                genomes.append([line, sequence])
                sequence = ""
            else:
                sequence = line + sequence
            if not len(lines) > 0:
                break
        return genomes
        
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
        description = (f"ORF starts at nucleotide {self.index+(self.frame%3)} and contains {self.length} nucleotides.\n"
            f"ORF Sequence: {self.sequence}\nThe Sequence is on the ")
        description += "antisense strand.\n" if (self.frame)>2 else "sense strand.\n"
        return description
        
    def is_longer(self, other_orf):
        return True if self.length >= other_orf.length else False
    
    def levenshtein_substitution_matrix(self, other_orf, insertion_penalty=1, deletion_penality=1, substituion_penalty=1):
        reference = self.sequence
        variance = other_orf.sequence
        substitution_matrix = [[0 for _ in range(len(variance)+1)] for _ in range(len(reference)+1)]
        
        for i in range(len(reference)+1):
            for j in range(len(variance)+1):
                if i==j==0:
                    substitution_matrix[i][j] = 0
                elif 1 <= j <= len(variance) and i==0:
                    substitution_matrix[i][j] = substitution_matrix[i][j-1] + insertion_penalty
                elif 1 <= i <= len(reference) and j==0:
                    substitution_matrix[i][j] = substitution_matrix[i-1][j] + deletion_penality
                else:
                    if reference[i-1] == variance[j-1]:
                        substitution_matrix[i][j] = substitution_matrix[i-1][j-1]
                    else:
                        substitution_matrix[i][j] = min([substitution_matrix[i][j-1]+insertion_penalty, substitution_matrix[i-1][j]+deletion_penality, substitution_matrix[i-1][j-1]+substituion_penalty])
        return substitution_matrix
    
    def find_snps(self, other_orf):
        sub_matrix = self.levenshtein_substitution_matrix(other_orf)
        snps = []
        i = j = 0
        print(sub_matrix)
        while i < len(self.sequence) or j < len(other_orf.sequence):
            if j == len(other_orf.sequence):
                snps.append([i, 'D', self.sequence[i]])
                i += 1
                continue
            if i == len(self.sequence):
                snps.append([i, 'I', other_orf.sequence[j]])
                j += 1
                continue
            the_lowest_distance = min(sub_matrix[i+1][j+1], sub_matrix[i][j+1], sub_matrix[i+1][j])
            if the_lowest_distance == sub_matrix[i+1][j+1]:
                if sub_matrix[i][j] != sub_matrix[i+1][j+1]:
                    snps.append([i+1, 'S', f"{self.sequence[i]}->{other_orf.sequence[j]}"])
                    # Find out if snp synonamous and potentially harmful by using two helper functions. is_sub_nonsym and is_sub_on_wanted_list
                i += 1
                j += 1
                continue
            if the_lowest_distance == sub_matrix[i][j+1]:
                snps.append([i+1, 'I', other_orf.sequence[j]])
                j += 1
                continue
            if the_lowest_distance == sub_matrix[i+1][j]:
                snps.append([i+1, 'D', self.sequence[i]])
                i += 1
                continue
        return snps
    
    def get_amino_acid(self):
        aa_sequence = ""
        for i in range(0, len(self.sequence), 3):
            codon = codon_table[self.sequence[i:i+3]]
            if codon != '*':
                aa_sequence += codon
        return aa_sequence
    
class AminoAcidChain:
    def __init__(self, acid_sequence="", orf=None):
        if orf:
            self.sequence = orf.get_amino_acid()
        else:
            self.sequence = acid_sequence
            
    def description(self):
        description = (f"This amino acid chain is {len(self.sequence)} acids long.\n"
                       f"Sequence: {self.sequence}")
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
        
    def give_sequence_length(self):
        return len(self.sequence)
    
    def get_longest_orf(self):
        longest_orf = ORF("", 0, 0, 0)
        for translation in range(6):
            if translation > 2:
                search_area = self.anti_sense
            else:
                search_area = self.sequence
            for i in range(translation%3, len(search_area), 3):
                head_codon = search_area[i:i+3]
                if len(head_codon) > 2 and codon_table[head_codon] == 'M':
                    for j in range(i+3, len(search_area) - 2, 3):
                        tail_codon = search_area[j:j+3]
                        if codon_table[tail_codon] == '*':
                            orf = ORF(search_area[i:j+3], i, j-i+3, translation)
                            if orf.is_longer(longest_orf):
                                longest_orf = orf
                            break
        return longest_orf
    
    def get_subsequence_as_ORF(self, start_index, length, translation):
        if translation > 2:
            search_area = self.anti_sense
        else:
            search_area = self.sequence
        return ORF(search_area[start_index:start_index+length], start_index, length, translation)