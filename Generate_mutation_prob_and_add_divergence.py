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
parser = argparse.ArgumentParser(description='Add divergence to reads in a FASTQ file based on a mutation matrix derived from an MSA.')
parser.add_argument('-i', '--input', required=True, help='Input compressed FASTA file (.fa.gz)')
parser.add_argument('-m', '--msa', required=True, help='Input MSA file (.fasta)')
parser.add_argument('-d', '--divergence', type=float, default=0.10, help='Divergence percentage (default: 0.10 for 10%%)')
parser.add_argument('-o', '--output', required=True, help='Output fasta file')

args = parser.parse_args()

"""
# Usage eg : 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Load_mutation_prob_and_add_divergence.py -i Simulation_Yersinia_pestis_SRR23219949_merged_100000_3DOperc_d.fa.gz\
 -m /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/aDNA_paper/dataset/Mappings/Yersinia_pestis_simulations/sibeliaz_alignment/Genome_to_MSA/output_sequences.fasta\
 -d 0.15
"""

fasta_file = args.input
divergence = args.divergence
MSA = args.msa

"""
MSA="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/aDNA_paper/dataset/Mappings/Yersinia_pestis_simulations/sibeliaz_alignment/Genome_to_MSA/output_sequences.fasta"
divergence = 0.15
fasta_file = "Simulation_Yersinia_pestis_SRR23219949_merged_100000_3DOperc_d.fa.gz"
"""


def calculate_mutation_matrix(fasta_file):
    """Calculate mutation probabilities from an aligned FASTA file."""
    # Initialize a dictionary to store mutation counts between each base
    mutation_counts = defaultdict(lambda: defaultdict(int))
    base_totals = defaultdict(int)  # To count total occurrences of each base
    # Read the aligned sequences from the FASTA file
    sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    # Ensure the sequences are of equal length (aligned)
    if len(set(len(seq) for seq in sequences)) != 1:
        raise ValueError("All sequences must be aligned and have the same length.")
    seq_length = len(sequences[0])
    # Compare each base at each position across all sequences
    for i in range(seq_length):
        column_bases = [seq[i] for seq in sequences]  # Bases at the i-th position across all sequences
        for j in range(len(column_bases)):
            base = column_bases[j]
            base_totals[base] += 1  # Count how often each base appears
            for k in range(len(column_bases)):
                if j != k and column_bases[j] != column_bases[k] and column_bases[k] != '-':  # Ignore gaps
                    mutation_counts[base][column_bases[k]] += 1
    # Convert counts to probabilities
    mutation_probs = {}
    for base, transitions in mutation_counts.items():
        total_mutations = sum(transitions.values())
        if total_mutations > 0:
            mutation_probs[base] = {mut_base: count / total_mutations for mut_base, count in transitions.items()}
        else:
            mutation_probs[base] = {}
    return mutation_probs

def mutate_base(base):
    if base in mutation_probs:
        mutations = list(mutation_probs[base].items())
        bases, probs = zip(*mutations)
        return random.choices(bases, probs)[0]
    return base  # Return the base unchanged if not in mutation_probs

# Load mutation probabilities from the MSA fil
#
##########

if MSA == "Yersinia":
    print("Mutation matrix already exists")
    mutation_probs = {
    "-": {
        "T": 0.26020863021743335,
        "G": 0.2403273442521625,
        "A": 0.2598767167558444,
        "C": 0.23951337007993664,
        "N": 7.39386946230839e-05
    },
    "G": {
        "C": 0.19266848488283408,
        "A": 0.6086973148593575,
        "T": 0.19862125689111707,
        "N": 1.29433666912656e-05
    },
    "C": {
        "G": 0.19213441153882485,
        "T": 0.6142538562586942,
        "A": 0.19359921585055162,
        "N": 1.2516351929307645e-05
    },
    "T": {
        "C": 0.6128790674177101,
        "G": 0.19762737178427722,
        "A": 0.18948458480469701,
        "N": 8.975993315616803e-06
    },
    "A": {
        "G": 0.6128069408903283,
        "T": 0.1917232599550518,
        "C": 0.19544808123135096,
        "N": 2.1717923268997607e-05
    },
    "N": {
        "A": 0.38461538461538464,
        "G": 0.23076923076923078,
        "C": 0.22377622377622378,
        "T": 0.16083916083916083
    }
    }

else:
	print(f"Reading MSA from {MSA} and calculating mutation probabilities...")
	mutation_probs = calculate_mutation_matrix(MSA) 
	print(f"Mutation probabilities: {json.dumps(mutation_probs, indent=4)}")


directory = os.path.dirname(fasta_file)
#output_filename = os.path.join(directory, f"{os.path.basename(fasta_file).replace('.fa.gz', f'_div-{int(divergence*100)}IDperc.fa.gz')}")
output_filename = args.output
print(f"Adding {divergence*100}% divergence to {fasta_file} based on the MSA...")
with gzip.open(fasta_file, "rt") as input_file, gzip.open(os.path.basename(output_filename), "wt") as output_file:
        for record in SeqIO.parse(input_file, "fasta"):  # Ensure you're reading the correct format
            mutated_seq = MutableSeq(record.seq)
            for i in range(len(mutated_seq)):
                if random.random() < divergence:
                    mutated_seq[i] = mutate_base(mutated_seq[i])  # Pass mutation_probs
            mutated_record = SeqRecord(mutated_seq, id=record.id, description=record.description)
            SeqIO.write(mutated_record, output_file, "fasta")
    
print(f"Mutated sequences written to {output_filename}")
