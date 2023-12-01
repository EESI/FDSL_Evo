import os
import pandas as pd
import pickle
import argparse


def combine_pickle_files_to_csv(folder_path, jobID, output_path):
    # Get a list of all pickle files in the folder
    pickle_files = [file for file in os.listdir(folder_path) if file.endswith(".pkl")]

    # Initialize an empty list to store the data frames
    data_frames = []

    # Loop through each pickle file, read it, and append its data frame to the list
    for pickle_file in pickle_files:
        file_path = os.path.join(folder_path, pickle_file)
        with open(file_path, "rb") as f:
            data_frame = pickle.load(f)
            data_frames.append(data_frame)

    # Combine all the data frames into one
    combined_df = pd.concat(data_frames, ignore_index=True)
    results = combined_df.drop_duplicates(subset=["SMILES"],keep='first')
    resultsSorted = results.sort_values(by = ['VINA'])

    # Write the combined data frame to a CSV file in the same folder
    output_csv_file = os.path.join(output_path, str(jobID)+"combined_data.csv")
    resultsSorted.to_csv(output_csv_file, index=False)


if __name__ == "__main__":
    argParse = argparse.ArgumentParser()

    argParse.add_argument("-j", "--job", help="SLURM Array Task ID used to generate result file name")
    argParse.add_argument("-p", "--pickledirectory", help="Path to directory where picklefiles are stored")
    argParse.add_argument("-o", "--outputdirectory", help="Path to directory where output CSV should be stored")

    args = argParse.parse_args()

    combine_pickle_files_to_csv(args.pickledirectory, args.job, args.outputdirectory)