import os
import subprocess
import argparse

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Run MAPS analysis on IBD data")

# Add arguments with a default value or not 
parser.add_argument("--input_folder", "-i", type=str, required=True, help="Path to the input folder")
parser.add_argument("--params_file", "-p", type=str, required=True, help="Path to the params file")
args = parser.parse_args()

# Save the prefix of the input file  
prefix = os.path.basename(os.path.normpath(args.input_folder))

## Update LD_LIBRARY_PATH environment variable
os.environ['LD_LIBRARY_PATH'] = "/home/eychenne@newcefe.newage.fr/Documents/BIOINFO/miniconda3/envs/MAPS/lib"

## Go to the folder 
os.chdir(args.input_folder)

# Run the command 
command = ['/home/eychenne@newcefe.newage.fr/Documents/BIOINFO/MAPS/src/runeems2', '--params', args.params_file]
process = subprocess.run(command, check=True)

print("MAPS run with success")

