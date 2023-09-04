import os
from openfold.utils.script_utils import parse_fasta


def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]


def load_fasta_sequences(fasta_dir):
    all_sequences = {}
    for fasta_file in list_files_with_extensions(fasta_dir, (".fasta", ".fa")):
        # Gather input sequences
        with open(os.path.join(fasta_dir, fasta_file), "r") as fp:
            data = fp.read()

        tags, seqs = parse_fasta(data)
        for tag, seq in zip(tags, seqs):
            all_sequences[tag] = seq
    return all_sequences
