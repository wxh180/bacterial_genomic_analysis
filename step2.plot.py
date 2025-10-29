import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

# Load your leader sequence CSV file
df = pd.read_csv("E.coli_tRNA_leader_sequences.csv")

# Analyze dinucleotides at positions -2 and -1 (last two bases of 5' leader)
dinucleotides = []

for seq in df["leader_seq"]:
    seq = seq.upper().replace('T', 'U')  # Convert DNA to RNA
    if len(seq) >= 2:
        dinuc = seq[-2:]  # last 2 bases are -2 and -1 positions
        if len(dinuc) == 2:
            dinucleotides.append(dinuc)

# Count dinucleotide frequencies
counts = Counter(dinucleotides)
total = sum(counts.values())

# Sort by frequency
sorted_counts = dict(sorted(counts.items(), key=lambda x: -x[1]))

# Plot donut chart
fig, ax = plt.subplots(figsize=(6, 6))
labels = list(sorted_counts.keys())
sizes = [v for v in sorted_counts.values()]
colors = plt.cm.tab20.colors[:len(labels)]

# Draw pie and white center for donut shape
wedges, texts = ax.pie(sizes, labels=labels, startangle=90, colors=colors)
centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig.gca().add_artist(centre_circle)

plt.title("Dinucleotide at positions −2 and −1\n(E. coli 5′ tRNA leaders)")
plt.tight_layout()
plt.show()
