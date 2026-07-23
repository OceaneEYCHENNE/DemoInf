import argparse  
import os

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to create params file for MAPS analysis")

# Add arguments with a default value or not 
parser.add_argument("--datapath", "-i", type=str,required=True, help="Path for data input")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output")
parser.add_argument("--output_filename", "-f", type=str, required=True, help="Name for output file")
parser.add_argument("--mcmcpath", "-m", type=str, required=True, help="Path for MCMC output")
parser.add_argument("--haplotypes", "-hap",type=int, default=200, help="Number of haplotypes sampled : 2 times number of individuals sampled")
parser.add_argument("--demes", "-d", type=int, default=100, help="Number of demes")
parser.add_argument("--genome_length", "-g", type=float, default=100.0, help="Value of genome size (in cM)")
parser.add_argument("--numMCMCIter", "-MCMC", type=int, default=3000000, help="Number MCMC iterations")
parser.add_argument("--numBurnIter", "-burn", type=int, default=1500000, help="Number burning iterations")
parser.add_argument("--numThinIter", "-thin", type=int, default=1500, help="Number of thinning iterations")
parser.add_argument("--lower_bound", "-lwr_b", type=str, default = '2', help="Value of lower bound (in cM)")
parser.add_argument("--upper_bound", "-upp_b", type=str, default = 'Inf', help="Value of upper bound (in cM)")
args = parser.parse_args()

## Create params file ## 
params = {
    "datapath": args.datapath,
    "mcmcpath": args.mcmcpath,
    "nIndiv": args.haplotypes,
    "nDemes": args.demes,
    "genomeSize": args.genome_length,
    "numMCMCIter": args.numMCMCIter,
    "numBurnIter": args.numBurnIter,
    "numThinIter": args.numThinIter,
    "lowerBound": args.lower_bound,
    "upperBound": args.upper_bound  
}

# Create output path for params file
output_folder = args.output_folder
output_filename = args.output_filename
output_path = os.path.join(output_folder, output_filename)

# Write file 
with open(output_path, "w") as f:
    for key, value in params.items():
        f.write(f"{key} = {value}\n")

print(f"Configuration file save here : {output_path}")
