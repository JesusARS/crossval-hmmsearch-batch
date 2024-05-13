import sys
import os
from Bio import SeqIO
import subprocess
import math


def get_word_size(identity_threshold):
    if 1 >= identity_threshold >= 0.7:
        return 5
    elif 0.7 > identity_threshold >= 0.6:
        return 4
    elif 0.6 > identity_threshold >= 0.5:
        return 3
    elif 0.5 > identity_threshold >= 0.4:
        return 2


def reduce_redundancy(file, job_name, redundancy):

    folder_name = f"{job_name}_{redundancy}"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    identity_threshold = redundancy/100
    output_name = f"{job_name}_{redundancy}.fasta"
    output_file = os.path.join(folder_name, output_name)

    word_size = get_word_size(identity_threshold)

    run_subprocess(f"cd-hit -i {file} -o {output_file} -c {identity_threshold} -n {word_size}")
    return folder_name, output_file


def extract_sequences(input_file, output_file, num_sequences):

    records_to_extract = []
    with open(input_file, 'r', encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records_to_extract.append(record)
            if len(records_to_extract) == num_sequences:
                break

    with open(output_file, 'w', encoding="utf-8") as output_handle:
        SeqIO.write(records_to_extract, output_handle, 'fasta')

    ids_to_extract = set(record.id for record in records_to_extract)

    with open(input_file, 'r', encoding="utf-8") as handle:
        records_to_remain = [record for record in SeqIO.parse(handle, 'fasta') if record.id not in ids_to_extract]

    with open(input_file, 'w', encoding="utf-8") as handle:
        SeqIO.write(records_to_remain, handle, 'fasta')


def k_fold_creation(folder_name, redundancy_file, k_folds):

    folds_folder = os.path.join(folder_name, 'Folds')
    if not os.path.exists(folds_folder):
        os.makedirs(folds_folder)

    num_sequences = count_sequences(redundancy_file)
    min_sequences_per_fold = math.floor(num_sequences/k_folds)
    num_folds_extra_sequence = num_sequences - min_sequences_per_fold*k_folds
    shuffle_file = os.path.join(folder_name, "redundancy_file_shuffle.fasta")

    run_subprocess(f"seqkit shuffle --rand-seed 1 {redundancy_file} -o {shuffle_file}")

    fold_count = 1

    for _ in range(1, num_folds_extra_sequence + 1):
        fold_file = os.path.join(folds_folder, f"fold_{fold_count}.fasta")
        sequences_number = min_sequences_per_fold + 1
        extract_sequences(shuffle_file, fold_file, sequences_number)
        fold_count += 1

    for _ in range(1, (k_folds - num_folds_extra_sequence) + 1):
        fold_file = os.path.join(folds_folder, f"fold_{fold_count}.fasta")
        extract_sequences(shuffle_file, fold_file, min_sequences_per_fold)
        fold_count += 1

    os.remove(shuffle_file)
 
    return folds_folder


def training_set(folds_folder, concat_file, iteration_name, fold):
    training_fold_path = os.path.join(folds_folder, f"fold_{fold}.fasta")
    result_file = f"result_{iteration_name}_{fold}.fasta"

    run_subprocess(f"touch {result_file}")
    run_subprocess(f"cat {training_fold_path} {concat_file} > {result_file}")

    os.remove(concat_file)

    return result_file


def training_and_testing_data_creation(folder_name, folds_folder, k_folds, negative_set_file):

    msa_models_folder = os.path.join(folder_name, 'MSA_and_Models')
    if not os.path.exists(msa_models_folder):
        os.makedirs(msa_models_folder)

    results_folder = os.path.join(folder_name, 'Results')
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    for i in range(1, k_folds + 1):
        iteration_name = f"msa_model_{i}"
        iteration_folder = os.path.join(msa_models_folder, iteration_name)
        if not os.path.exists(iteration_folder):
            os.makedirs(iteration_folder)

        testing_set_file_path = os.path.join(folds_folder, f"fold_{i}.fasta")
        testing_set_copy_file_path = os.path.join(iteration_folder, f"testing_set_{i}.fasta")
        run_subprocess(f"cp {testing_set_file_path} {testing_set_copy_file_path}")
        
        concat_file = "empty.fasta"
        run_subprocess(f"touch {concat_file}")

        for fold in range(1, k_folds + 1):
            if fold != i:
                concat_file = training_set(folds_folder, concat_file, iteration_name, fold)

        training_set_file_path = os.path.join(iteration_folder, f"training_set_{i}.fasta")
        run_subprocess(f"touch {training_set_file_path}")
        run_subprocess(f"cp {concat_file} {training_set_file_path}")

        os.remove(concat_file)

        msa_path = os.path.join(iteration_folder, f"{iteration_name}_msa.sto")
        run_subprocess(f"./clustalo -i {training_set_file_path} -o {msa_path} --outfmt=st")

        hmm_path = os.path.join(iteration_folder, f"{iteration_name}_hmm.hmm")
        run_subprocess(f"hmmbuild {hmm_path} {msa_path}")

        result_path = os.path.join(results_folder, f"{iteration_name}_hmm_result")
        run_subprocess(f"cat {testing_set_copy_file_path} {negative_set_file} > full_testing_data.fasta")
        run_subprocess(f"hmmsearch --max --tblout {result_path} {hmm_path} full_testing_data.fasta")

        os.remove("full_testing_data.fasta")


def msa_creation(iteration_folder, training_set_file_path):
    msa_path = os.path.join(iteration_folder, 'msa.sto')
    run_subprocess(f"./clustalo -i {training_set_file_path} -o {msa_path} --outfmt=st")

def count_sequences(file):

    num_sequences = 0

    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            num_sequences += 1

    return num_sequences


def run_subprocess(command):
    subprocess.run(command, shell=True, capture_output=True, text=True, check=False)


def main():

    job_name = sys.argv[1]
    training_file = sys.argv[2]
    k_folds = sys.argv[3]
    redundancy = sys.argv[4]
    negative_set_file = sys.argv[5]
    redundancy_str = redundancy.split(',')
    redundancy_list = [int(val) for val in redundancy_str]

    for redundancy in redundancy_list:
        folder_name, redundancy_file = reduce_redundancy(training_file, job_name, redundancy)
        folds_folder = k_fold_creation(folder_name, redundancy_file, int(k_folds))
        training_and_testing_data_creation(folder_name, folds_folder, int(k_folds), negative_set_file)


if __name__ == "__main__":
    main()
