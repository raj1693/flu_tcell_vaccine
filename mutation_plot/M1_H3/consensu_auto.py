import os
import subprocess
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

# Define base folders
base_folder = "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/M1_H3/"  # Replace with your actual clade folder path
output_folder = "/Users/raj/Desktop/Codes_bioinfro/R scripts/mutation_plot/M1_H3/python_output"

# Ensure output folder exists
os.makedirs(output_folder, exist_ok=True)

# Path to MAFFT executable (ensure it's installed and accessible)
mafft_exe = "mafft"  # If not in system PATH, provide full path like "/usr/bin/mafft"

# Process each clade separately
for clade in os.listdir(base_folder):
    clade_path = os.path.join(base_folder, clade)
    
    if os.path.isdir(clade_path):  # Ensure it's a directory
        merged_fasta = os.path.join(output_folder, f"{clade}_merged.fasta")
        aligned_fasta = os.path.join(output_folder, f"{clade}_aligned.fasta")
        consensus_output = os.path.join(output_folder, f"{clade}_consensus.txt")

        # Step 1: Merge all FASTA files within the clade
        with open(merged_fasta, "w") as outfile:
            for file in os.listdir(clade_path):
                if file.endswith(".fasta") or file.endswith(".fa"):
                    with open(os.path.join(clade_path, file), "r") as infile:
                        outfile.write(infile.read())

        # Step 2: Align sequences using MAFFT
        with open(aligned_fasta, "w") as aligned_out:
            subprocess.run([mafft_exe, "--auto", merged_fasta], stdout=aligned_out, check=True)

        # Step 3: Compute Consensus Sequence
        alignment = AlignIO.read(aligned_fasta, "fasta")
        summary = SummaryInfo(alignment)
        consensus_seq = summary.dumb_consensus()

        # Step 4: Save consensus sequence
        with open(consensus_output, "w") as out_file:
            out_file.write(f">Consensus_{clade}\n{consensus_seq}\n")

        print(f"Consensus for {clade} saved at {consensus_output}")

