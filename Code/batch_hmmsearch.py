"""Module that executes hmmsearch in batches."""
import argparse
import csv
import os
import subprocess
from Bio import SeqIO
from HmmDataProcessor import HmmDataProcessor
import pandas as pd


def run_subprocess(command):
    """
    Run a shell command as a subprocess.

    Args:
        command (str): The shell command to execute.
    """

    subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True,
        check=False
    )


def quantificate_hmm_data(results_folder, hmm_csv_path, fasta_label):

    df = pd.read_csv(hmm_csv_path)

    # Append quantification to the CSV file
    output_file = os.path.join(results_folder, "quantification.csv")

    with open(output_file, "a", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "label", "quantity",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header if the file is empty
        if csvfile.tell() == 0:
            writer.writeheader()

        # Write the row of results
        writer.writerow({
            "label": fasta_label,
            "quantity": df['target_name'].nunique(),
        })


def perform_hmm_search(hmm_config, fasta_label, fasta_path):

    results_subfolder = os.path.join(
        hmm_config["results_folder"], fasta_label.replace(' ', '_')
    )
    os.makedirs(results_subfolder, exist_ok=True)

    redundancy_path = os.path.join(
        results_subfolder, f"{fasta_path}_100"
    )
    try:
        run_subprocess(
            f"cd-hit -i {fasta_path} -o {redundancy_path} -c 1 -n 5"
        )
    except Exception as e:
        raise RuntimeError(f"Error running CD-HIT: {e}") from e

    results_path = os.path.join(
        results_subfolder, f"{fasta_label.replace(' ', '_')}_hmm_result"
    )
    hmm_csv_path = os.path.join(
        results_subfolder, f"{fasta_label.replace(' ', '_')}_hmm_table.csv"
    )
    proteins_path = os.path.join(
        results_subfolder, "proteins.fasta"
    )

    run_subprocess(
        f"hmmsearch --max --seed 1 -E {hmm_config['threshold']} --tblout {results_path} {hmm_config['hmm_path']} {redundancy_path}"
    )

    HmmDataProcessor().parse_results_data(results_path, hmm_csv_path)
    quantificate_hmm_data(
        hmm_config["results_folder"], hmm_csv_path, fasta_label
    )

    df = pd.read_csv(hmm_csv_path)
    
    with open(proteins_path, "w", encoding="utf-8") as output_handle:
        for protein_id in df["target_name"]:
            for record in SeqIO.parse(fasta_path, "fasta"):
                if record.id == protein_id:
                    SeqIO.write(record, output_handle, "fasta")

    run_subprocess(
        f"cat {proteins_path} >> {hmm_config['all_sequences_path']}"
    )


def main():

    description = """
    Script that runs hmmsearch on a batch of FASTA files.

    This script takes an HMM file and runs hmmsearch on a set of FASTA files
    specified in a CSV file. The CSV file must have a first column indicating
    the file label, and a second column indicating its path.

    A 'results' folder will be created to store the result table for each FASTA
    file, and a 'quantification.csv' file will be generated to indicate the
    number of proteins identified in each file, based on a specified E-value
    threshold.


    System Requirements:
    - Linux operating system
    - External software dependencies:
        + hmmer (http://hmmer.org/)

    Note: Install hmmer and ensure it is accessible in your system PATH.
    """

    epilog = """
    Example usage:
    --------------
    python3 batch_hmmsearch.py GAF HMM.hmm FASTA_files.csv 1e-5
    """

    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "job_name",
        type=str,
        help="Job name"
    )
    parser.add_argument(
        "hmm_path",
        type=argparse.FileType("r"),
        help="HMM path"
    )
    parser.add_argument(
        "FASTA_files",
        type=argparse.FileType("r"),
        help="CSV file with the path of the FASTA files"
    )
    parser.add_argument(
        "E_value",
        type=float,
        help="E-value threshold"
    )

    args = parser.parse_args()

    results_folder = args.job_name
    os.makedirs(results_folder, exist_ok=True)
    df = pd.read_csv(args.FASTA_files.name)
    all_sequences_path = os.path.join(results_folder, "sequences.fasta")

    run_subprocess(
        f"touch {all_sequences_path}"
    )

    hmm_config = {
        "hmm_path": args.hmm_path.name,
        "all_sequences_path": all_sequences_path,
        "threshold": args.E_value,
        "results_folder": results_folder,
    }

    for _, row in df.iterrows():

        fasta_label = row.iloc[0]
        fasta_path = row.iloc[1]
        perform_hmm_search(hmm_config, fasta_label, fasta_path)


if __name__ == "__main__":
    main()
