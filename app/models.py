nucleotide_complements = {
    'A': 'T', 'T': 'A',
    
    'G': 'C', 'C': 'G',
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
        
        
    
    
class Sequence:
    """
    This is a model of DNA and mRNA sequences. This model wraps up functions
    useful for processing these sequences.
    """
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        
    def give_sequence_length(self):
        return len(self.sequence)
        
    def create_antisense(self):
        anti_sense_s = ""
        for char in reversed(self.sequence):
            anti_sense_s += nucleotide_complements[char]
        return(anti_sense_s)
        