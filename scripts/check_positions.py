"""
For training purposes. To find out which aminoacids in the tested protein are occupying selected positions.
"""

from Bio import SeqIO

# Načteme sekvenci
record = next(SeqIO.parse('../data/example_protein.fasta', 'fasta'))
seq = str(record.seq)

# Pozice z našeho CSV (nezapomeň, že v biologii se počítá od 1)
positions = [7, 22, 45, 64, 103]

print("Checking amino acids at specified positions:")
print("Position | Amino Acid")
print("-" * 25)

for pos in positions:
    if pos <= len(seq):
        aa = seq[pos - 1]  # -1 protože Python indexuje od 0
        print(f"{pos:8} | {aa}")
    else:
        print(f"{pos:8} | OUT OF RANGE (sequence has {len(seq)} amino acids)")