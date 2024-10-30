from app.models import FileReader, Sequence

fr = FileReader()

text_sequences = fr.extract_sequences("gene.fna")

sequence_object = []
for sq in text_sequences:
    sequence_object.append(Sequence(sq[0], sq[1]))

for sq in sequence_object:
    print(sq.create_antisense())