"""Module that executes k-fold cross-validation on protein sequence data"""
import argparse
import csv
import os
import subprocess
from Bio import SeqIO
from HmmDataProcessor import HmmDataProcessor
import pandas as pd


def get_word_size(identity_threshold):
    """
    Determine the word size based on the identity threshold.

    Args:
        identity_threshold (float): The identity threshold, which should be
            less than or equal to 1.0.

    Returns:
        int: The word size corresponding to the given identity threshold.
    """
    if 1 >= identity_threshold >= 0.7:
        return 5
    if 0.7 > identity_threshold >= 0.6:
        return 4
    if 0.6 > identity_threshold >= 0.5:
        return 3
    return 2


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


def reduce_redundancy(training_data_file, job_name, redundancy):
    """
    Reduce redundancy in a sequence file using CD-HIT and save the output to a
    new file.

    Args:
        training_data_file (str): Path to the input sequence file with the
            training data.
        job_name (str): The job name used to create the output folder and
            file name.
        redundancy (int): The redundancy level as a percentage (0 to 100).

    Returns:
        dict: A dictionary containing the folder name, the path to the
            redundancy file, and the redundancy number.
    """

    folder_name = f"{job_name}_{redundancy}"
    os.makedirs(folder_name, exist_ok=True)
  
    identity_threshold = redundancy / 100
    output_name = f"{job_name}_{redundancy}.fasta"
    redundancy_file = os.path.join(folder_name, output_name)

    word_size = get_word_size(identity_threshold)

    try:
        run_subprocess(
            f"cd-hit -i {training_data_file} -o {redundancy_file} -c {identity_threshold} -n {word_size}"
        )
    except Exception as e:
        raise RuntimeError(f"Error running CD-HIT: {e}") from e

    return {
        "folder_name": folder_name,
        "redundancy_file": redundancy_file,
        "redundancy": redundancy,
    }


def extract_sequences(input_file, output_file, num_sequences):
    """
    Extract a specified number of sequences from a FASTA file and update
        the input file.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file where extracted
            sequences will be saved.
        num_sequences (int): Number of sequences to extract.

    Returns:
        None
    """
    records_to_extract = []
    records_to_remain = []
    
    with open(input_file, 'r', encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if len(records_to_extract) < num_sequences:
                records_to_extract.append(record)
            else:
                records_to_remain.append(record)
    
    with open(output_file, 'w', encoding="utf-8") as output_handle:
        SeqIO.write(records_to_extract, output_handle, 'fasta')

    with open(input_file, 'w', encoding="utf-8") as handle:
        SeqIO.write(records_to_remain, handle, 'fasta')


def k_fold_creation(redundancy_config, k_folds):
    """
    Create k folds for cross-validation from a redundancy-filtered sequence
    file.

    Args:
        redundancy_config (dict): A dictionary with the following items:
            folder_name (str): Directory where the folds will be created.
            redundancy_file (str): Path to the redundancy-filtered sequence
                file.

        k_folds (int): Number of folds to create.

    Returns:
        str: Path to the directory containing the generated folds.
    """

    folder_name = redundancy_config["folder_name"]
    redundancy_file = redundancy_config["redundancy_file"]

    folds_folder = os.path.join(folder_name, "folds")
    os.makedirs(folds_folder, exist_ok=True)

    num_sequences = count_sequences(redundancy_file)
    min_sequences_per_fold = num_sequences // k_folds
    num_folds_extra_sequence = num_sequences % k_folds
    shuffle_file = os.path.join(folder_name, "redundancy_file_shuffle.fasta")

    try:
        run_subprocess(
            f"seqkit shuffle --rand-seed 1 {redundancy_file} -o {shuffle_file}"
        )
    except Exception as e:
        raise RuntimeError(f"Failed to shuffle sequences: {e}") from e

    for iteration in range(1, k_folds + 1):
        fold_file = os.path.join(folds_folder, f"fold_{iteration}.fasta")
        if iteration <= num_folds_extra_sequence:
            sequences_number = min_sequences_per_fold + 1
        else:
            sequences_number = min_sequences_per_fold
        extract_sequences(shuffle_file, fold_file, sequences_number)

    try:
        os.remove(shuffle_file)
    except Exception as e:
        raise RuntimeError(
            f"Failed to remove temporary shuffle file: {e}"
        ) from e

    return folds_folder


def create_empty_file(file_path):
    """Create an empty file at the specified path."""
    run_subprocess(f"touch {file_path}")


def copy_file(source, destination):
    """Copy a file from source to destination."""
    run_subprocess(f"cp {source} {destination}")


def concatenate_files(input_file_1, input_file_2, output_file):
    """Concatenate two files into one."""
    run_subprocess(f"cat {input_file_1} {input_file_2} > {output_file}")


def update_training_set(iteration, folds_folder, concat_file, fold):
    """
    Create the training set by concatenating the fold file with the current
    concatenated file.
   
    Args:
        folds_folder (str): Path to the folder containing the fold files.
        concat_file (str): Path to the current concatenated file.
        iteration (str): Name of the current iteration.
        fold (int): Fold number to process.

    Returns:
        str: Path to the resulting training set file.
    """
    training_fold_path = os.path.join(folds_folder, f"fold_{fold}.fasta")
    result_file = f"result_{iteration}_{fold}.fasta"

    create_empty_file(result_file)
    concatenate_files(training_fold_path, concat_file, result_file)

    os.remove(concat_file)

    return result_file


def training_and_testing_data_creation(redundancy_config, iteration, iteration_folder):
    """
    Create the training and testing datasets for the k-fold cross-validation
    iteration.

    Args:
        redundancy_config (dict): A dictionary with the following items:
            folds_folder (str): Path to the folder containing the fold files.
            k_folds (int): Number of folds for cross-validation.

        iteration (int): The current iteration number.
        iteration_folder (str): Path to the folder where the iteration results
            will be stored.

    Returns:
        str: Path to the training and testing set files for this iteration.
    """

    folds_folder = redundancy_config["folds_folder"]

    # Paths for testing set files
    testing_set_file_path = os.path.join(
        folds_folder, f"fold_{iteration}.fasta"
    )
    testing_set_copy_file_path = os.path.join(
        iteration_folder, f"testing_set_{iteration}.fasta"
    )
    
    # Copy the fold that corresponds to the current iteration to be the
    # testing set
    copy_file(testing_set_file_path, testing_set_copy_file_path)
    
    # Create an empty file to start concatenating the training set
    concat_file = "empty.fasta"
    create_empty_file(concat_file)

    # Iterate over all folds to create the training set by concatenating all
    # folds except the current testing fold
    for fold in range(1, redundancy_config["k_folds"] + 1):
        if fold != iteration:
            concat_file = update_training_set(
                iteration, folds_folder, concat_file, fold
            )

    # Save the final training set after all folds have been concatenated
    training_set_file_path = os.path.join(
        iteration_folder, f"training_set_{iteration}.fasta"
    )
    copy_file(concat_file, training_set_file_path)

    # Remove the temporary concatenation file
    os.remove(concat_file)

    return training_set_file_path, testing_set_copy_file_path


def build_msa_and_hmm(training_set, iteration_folder, iteration_name):
    """
    Generate a Multiple Sequence Alignment (MSA) and build a Hidden Markov
    Model (HMM).

    Args:
        training_set (str): Path to the input training sequences.
        iteration_folder (str): Directory to save the MSA and HMM files.
        iteration_name (str): Identifier for the iteration.

    Returns:
        str: Path to the generated HMM file.
    """
    msa_path = os.path.join(iteration_folder, f"{iteration_name}_msa.fasta")
    run_subprocess(f"clustalo -i {training_set} -o {msa_path} --outfmt=fasta")

    hmm_path = os.path.join(iteration_folder, f"{iteration_name}_hmm.hmm")
    run_subprocess(f"hmmbuild --seed 1 {hmm_path} {msa_path}")
  
    return hmm_path


def perform_hmm_search(hmm_config, negative_set, results_folder):
    """
    Run an HMM search on the combined testing and negative sets.

    Args:
        hmm_config (dict): A dictionary with the following items:
            hmm_path (str): Path to the HMM file.
            iteration (int): The current iteration number for the process.
            iteration_name (str): Identifier for the iteration.
            positive_control_file (str): Path to the FASTA file containing
                positive control protein sequences.

        negative_set (str): Path to the negative sequences.
        results_folder (str): Directory to save results.

    Returns:
        tuple: Path to the CSV results file.
    """
    results_path = os.path.join(
        results_folder, f"{hmm_config['iteration_name']}_hmm_result"
    )
    testing_data_path = "full_testing_data.fasta"

    concatenate_files(
        hmm_config["positive_control_file"], negative_set, testing_data_path
    )
    run_subprocess(
        f"hmmsearch --max --seed 1 --tblout {results_path} {hmm_config['hmm_path']} {testing_data_path}"
    )

    hmm_csv_path = os.path.join(
        results_folder, f"{hmm_config['iteration']}_hmm_table.csv"
    )
    HmmDataProcessor().parse_results_data(results_path, hmm_csv_path)

    os.remove(testing_data_path)

    return hmm_csv_path


def extract_statistical_data(hmm_config, target_file, total_negative_proteins):
    """
    Processes  HMMER search data to compute true positives (TP), false
    positives (FP), false negatives (FN), and true negatives (TN) based on
    different E-value thresholds. Results are appended to a CSV file.

    Args:
        hmm_config (dict): A dictionary with the following items:
                iteration (int): The current iteration number for the process.
                positive_control_file (str): Path to the FASTA file containing
                    positive control protein sequences.
                redundancy (str): Redundancy level used in the process.

        target_file (str): Path to the CSV file containing target protein data.
        total_negative_proteins (int): Total number of negative control
            proteins.

    Returns:
        None
    """

    # Define E-value thresholds from 10^0 to 10^-10
    e_value_thresholds = [10**-x for x in range(11)]

    for e_value_threshold in e_value_thresholds:
        # Read target protein data
        target_data = pd.read_csv(target_file)

        # Filter proteins with E-value below the threshold
        filtered_targets = target_data[
            target_data["E-value_fs"] < e_value_threshold
        ]

        # Read positive control protein IDs
        positive_control_ids = set()
        with open(hmm_config["positive_control_file"], "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                positive_control_ids.add(record.id)

        # Count true positives (TP) and false positives (FP)
        tp = filtered_targets["target_name"].isin(positive_control_ids).sum()
        fp = len(filtered_targets) - tp

        # Count false negatives (FN)
        total_positive_proteins = len(positive_control_ids)
        fn = total_positive_proteins - tp

        # Count true negatives (TN)
        tn = total_negative_proteins - fp

        # Append results to the CSV file
        output_file = "results.csv"
        with open(output_file, "a", newline="", encoding="utf-8") as csvfile:
            fieldnames = [
                "redundancy", "iteration",
                "E-value", "TP", 
                "TN", "FP", "FN"
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            # Write header if the file is empty
            if csvfile.tell() == 0:
                writer.writeheader()

            # Write the row of results
            writer.writerow({
                "redundancy": hmm_config["redundancy"],
                "iteration": hmm_config["iteration"],
                "E-value": e_value_threshold,
                "TP": tp,
                "TN": tn,
                "FP": fp,
                "FN": fn
            })


def models_creation(redundancy_config, negative_set_file):
    """
    Create models for each fold in k-fold cross-validation.

    Args:
        redundancy_config (dict): A dictionary with the following items:
            folder_name (str): The main directory where results and models
            will be stored.
            folds_folder (str): Directory containing the fold files.
            k_folds (int): Number of folds for cross-validation.
            redundancy (int): The redundancy level for the dataset.

        negative_set_file (str): Path to the negative testing data in fasta
            format.

    Returns:
        None
    """

    folder_name = redundancy_config["folder_name"]

    msa_models_folder = os.path.join(folder_name, 'MSA_and_models')
    os.makedirs(msa_models_folder, exist_ok=True)

    results_folder = os.path.join(folder_name, 'results')
    os.makedirs(results_folder, exist_ok=True)

    number_negative_sequences = count_sequences(negative_set_file)

    for iteration in range(1, redundancy_config["k_folds"] + 1):
        iteration_name = f"msa_model_{iteration}"
        iteration_folder = os.path.join(msa_models_folder, iteration_name)
        os.makedirs(iteration_folder, exist_ok=True)

        training_set, testing_set = training_and_testing_data_creation(
            redundancy_config, iteration, iteration_folder,
        )

        # Build MSA and HMM
        hmm_path = build_msa_and_hmm(
            training_set, iteration_folder, iteration_name
        )

        hmm_config = {
            "hmm_path": hmm_path,
            "iteration": iteration,
            "iteration_name": iteration_name,
            "positive_control_file": testing_set,
            'redundancy': redundancy_config["redundancy"],
        }
        # Perform HMM search and extract statistic data
        hmm_csv_path = perform_hmm_search(
            hmm_config, negative_set_file, results_folder
        )

        extract_statistical_data(
            hmm_config, hmm_csv_path,
            number_negative_sequences,
        )


def count_sequences(file):
    """
    Count the number of sequences in a FASTA file.

    Args:
        file (str): Path to the FASTA file.

    Returns:
        int: The number of sequences in the file.
    """
    with open(file, "r", encoding="utf-8") as handle:
        count = sum(1 for _ in SeqIO.parse(handle, "fasta"))
    return count


def main():
    """
    Execute k-fold cross-validation on protein sequence data with varying
    redundancy levels.
    """

    description = """
    Script that executes a k-fold cross-validation for protein sequence data.

    This script takes a training data file in FASTA format and performs k-fold
    cross-validation with different levels of redundancy. The redundancy levels
    are provided as a comma-separated list and will be evaluated as a
    hyperparameter. A negative set file is also required for testing purposes.

    System Requirements:
    - Linux operating system
    - External software dependencies:
        + hmmer (http://hmmer.org/)
        + clustalo (http://www.clustal.org/omega/)
        + seqkit (https://bioinf.shenwei.me/seqkit/)
    - Python packages:
        + Biopython (https://biopython.org/)

    Note: Install hmmer, clustalo and seqkit separately and ensure they are
    accessible in your system PATH.
    Ensure Biopython is installed in your Python environment.
    """

    epilog = """
    Example usage:
    --------------
    python3 k_fold_cross_validation.py GAF_domain -ts GAF_training_data.fasta -ns negative_data.fasta -k 5 -r 100,90,80,70,60
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
        "-ts", "--training_data_file",
        type=argparse.FileType("r"),
        required=True,
        help="Training data in FASTA format"
    )
    parser.add_argument(
        "-ns", "--negative_set_file",
        type=argparse.FileType("r"),
        help="Negative testing data in fasta format"
    )
    parser.add_argument(
        "-k", "--k_folds",
        type=int,
        default="5",
        help="Number of k fold cross validation"
    )
    parser.add_argument(
        "-r", "--redundancy",
        type=str,
        default="100",
        help="""Redudancy levels to be evaluated as a hyperparameter. Set a
                sequence of numbers from 100 to 40 separated with comma.
                Example: 100,90,80,70"""
    )

    args = parser.parse_args()

    try:
        redundancy_list = [int(val) for val in args.redundancy.split(',')]
    except ValueError as e:
        parser.error(f"Invalid redundancy levels provided: {e}")

    for redundancy in redundancy_list:

        redundancy_config = reduce_redundancy(
            args.training_data_file.name, args.job_name, redundancy
        )
        folds_folder = k_fold_creation(redundancy_config, args.k_folds)
        redundancy_config.update({
            "folds_folder": folds_folder,
            "k_folds": args.k_folds,
        })
        models_creation(
            redundancy_config, args.negative_set_file.name,
        )


if __name__ == "__main__":
    main()
