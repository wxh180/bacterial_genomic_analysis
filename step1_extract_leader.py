from Bio import SeqIO
import pandas as pd

GENOME_FASTA = 'GCF_000005845.2_ASM584v2_genomic.fna'
GFF_FILE = 'GCF_000005845.2_ASM584v2_genomic.gff'
LEADER_LENGTH = 2 
SPECIES = "E.coli"

genome = SeqIO.to_dict(SeqIO.parse(GENOME_FASTA, "fasta"))

trna_records = []

with open(GFF_FILE) as gff:
    for line in gff:
        if line.startswith("#") or "\ttRNA\t" not in line:
            continue
        parts = line.strip().split("\t")
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attributes = parts[8]
        gene = ""
        if "Name=" in attributes:
            gene = attributes.split("Name=")[1].split(";")[0]

        if chrom not in genome:
            continue
        seq = genome[chrom].seq
        if strand == "+":
            leader_start = max(0, start - LEADER_LENGTH - 1)
            leader_seq = seq[leader_start:start - 1]
        else:
            leader_end = min(len(seq), end + LEADER_LENGTH)
            leader_seq = seq[end:leader_end].reverse_complement()

        trna_records.append({
            'species': SPECIES,
            'gene': gene,
            'strand': strand,
            'leader_seq': str(leader_seq)
        })

df = pd.DataFrame(trna_records)
df.to_csv(f"{SPECIES}_tRNA_leader_sequences.csv", index=False)
print(f"Saved to {SPECIES}_tRNA_leader_sequences.csv")
