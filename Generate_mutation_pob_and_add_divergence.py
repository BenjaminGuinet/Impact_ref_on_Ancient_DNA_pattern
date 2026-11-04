import argparse
import gzip
import random
import os
from Bio import SeqIO
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq


print("\n")

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Add divergence to a FASTA file based on a mutation matrix derived from an MSA.')
parser.add_argument('-i', '--input',      required=True,              help='Input FASTA file (plain or .gz)')
parser.add_argument('-m', '--msa',        required=True,              help='Input MSA file (plain or .gz)')
parser.add_argument('-d', '--divergence', type=float, default=0.10,   help='Divergence fraction (default: 0.10 for 10%%)')
parser.add_argument('-o', '--output',     required=True,              help='Output FASTA file (plain or .gz)')

args = parser.parse_args()

fasta_file = args.input
divergence = args.divergence
MSA        = args.msa
output_filename = args.output


# ------------------------------------------------------------------
# Helper: open plain or gzipped files transparently
# ------------------------------------------------------------------
def smart_open(path, mode="r"):
    """Open a file in text mode regardless of whether it is gzip-compressed."""
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t" if "t" not in mode else mode)
    # plain open – strip any accidental 't' suffix so it works for both cases
    plain_mode = mode.replace("t", "") if not path.endswith(".gz") else mode
    return open(path, plain_mode)


# ------------------------------------------------------------------
# Mutation-matrix calculation from MSA
# ------------------------------------------------------------------
def calculate_mutation_matrix(msa_file):
    """Calculate mutation probabilities from an aligned FASTA file."""
    mutation_counts = defaultdict(lambda: defaultdict(int))
    base_totals     = defaultdict(int)

    with smart_open(msa_file, "r") as fh:
        sequences = [str(record.seq) for record in SeqIO.parse(fh, "fasta")]

    if len(set(len(seq) for seq in sequences)) != 1:
        raise ValueError("All sequences in the MSA must be aligned and have the same length.")

    seq_length = len(sequences[0])
    for i in range(seq_length):
        column_bases = [seq[i] for seq in sequences]
        for j, base in enumerate(column_bases):
            base_totals[base] += 1
            for k, other in enumerate(column_bases):
                if j != k and base != other and other != '-':
                    mutation_counts[base][other] += 1

    mutation_probs = {}
    for base, transitions in mutation_counts.items():
        total = sum(transitions.values())
        if total > 0:
            mutation_probs[base] = {mb: cnt / total for mb, cnt in transitions.items()}
        else:
            mutation_probs[base] = {}
    return mutation_probs


# ------------------------------------------------------------------
# Mutate a single base according to pre-computed probabilities
# ------------------------------------------------------------------
def mutate_base(base, mutation_probs):
    if base in mutation_probs and mutation_probs[base]:
        bases, probs = zip(*mutation_probs[base].items())
        return random.choices(bases, probs)[0]
    return base


# ------------------------------------------------------------------
# Load / hard-code mutation probabilities
# ------------------------------------------------------------------

print(f"Reading MSA from {MSA} and calculating mutation probabilities...")
mutation_probs = calculate_mutation_matrix(MSA)
print(f"Mutation probabilities:\n{json.dumps(mutation_probs, indent=4)}")


# ------------------------------------------------------------------
# Main: apply divergence to every sequence in the input FASTA
# ------------------------------------------------------------------
print(f"Adding {divergence * 100:.2f}% divergence to {fasta_file}...")

n_records  = 0
n_mutated  = 0

with smart_open(fasta_file, "r") as input_fh, \
     smart_open(output_filename, "w") as output_fh:

    for record in SeqIO.parse(input_fh, "fasta"):
        n_records += 1
        mutated_seq = MutableSeq(str(record.seq))

        for i in range(len(mutated_seq)):
            if random.random() < divergence:
                new_base = mutate_base(mutated_seq[i], mutation_probs)
                if new_base != mutated_seq[i]:
                    n_mutated += 1
                mutated_seq[i] = new_base

        mutated_record = SeqRecord(mutated_seq,
                                   id=record.id,
                                   description=record.description)
        SeqIO.write(mutated_record, output_fh, "fasta")

print(f"Done. Processed {n_records} sequence(s), introduced {n_mutated} substitution(s).")
print(f"Output written to: {output_filename}")