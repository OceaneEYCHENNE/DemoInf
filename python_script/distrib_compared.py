import argparse  
import pandas as pd
import matplotlib 
from matplotlib import pyplot as plt

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to compare length distribution of IBD segments from two csv file")

# Add arguments with a default value or not 
parser.add_argument("--input_file1", "-i1", type=str,required=True, help="Path for IBD input file 1")
parser.add_argument("--input_file2", "-i2", type=str,required=True, help="Path for IBD input file 2")
parser.add_argument("--output_folder","-o", type=str, required=True, help="Path for output folder")
parser.add_argument("--chr_length", "L", type=int, default = 1e8, help="Length of simulated chromosome")
parser.add_argument("--min_length", "-l_min", type=int, default=None, help="Minimum length of represented segments")
parser.add_argument("--max_length", "-l_max", type=int, default=None, help="Maximum length of represented segments")
parser.add_argument("--min_time" ,"-t_min", type=int, default=None, help="Minimum age of represented segments")
parser.add_argument("--max_time", "-t_max", type=int, default=None, help="Maximum age of represented segments")
args = parser.parse_args()

# Save the prefix of the input file 1 
input_path1 = Path(args.input_file1)
filename1 = input_path1.name.split(".csv")[0]
prefix1 = filename1.split("IBD_coord_")[1]

# Save the prefix of the input file 2 
input_path2 = Path(args.input_file2)
filename2 = input_path2.name.split(".csv")[0]
prefix2 = filename2.split("IBD_coord_")[1]

### Import csv ###
dfIBD1 = pd.read_csv(args.input_file1)
dfIBD2 = pd.read_csv(args.input_file2)


### Filter IBD segments by minimal and/or maximal ages, if specified ### 
if args.min_time is not None:
	dfIBD1 = dfIBD1[dfIBD1.Time > args.min_time]
	dfIBD2 = dfIBD2[dfIBD2.Time > args.min_time]
	
if args.max_time is not None:
	dfIBD1 = dfIBD1[dfIBD1.Time < args.max_time]
	dfIBD2 = dfIBD2[dfIBD2.Time < args.max_time]

### Return the number of IBD segments recorded for each file ###
num_ibd_1 = len(dfIBD1.Pair) 
num_ibd_2 = len(dfIBD2.Pair) 
print(f"Loaded two csv files containing {num_ibd_1} IBD segments and {num_ibd_2} respectively")
ts_sim

### Retrieve variable values ###
L = args.chr_length

	
### Create comparison histogram of IBD segments length ### 
plt.figure()
plt.hist(dfIBD1["Length"], range = (0, L), bins=60, color = '#66CC99', edgecolor = 'white', alpha = 0.7, label = "data1")
plt.hist(dfIBD2["Length"], range = (0, L), bins=60, color = '#009999', edgecolor = 'white', alpha = 0.7, label = "data2")
plt.xlabel('Length')
plt.ylabel('Number')
plt.title(f'Length distribution of simulated IBD segments', fontsize = 16, loc = 'center')
plt.yscale('log')
plt.xlim(args.min_length, args.max_length)
plt.legend(loc = 'upper right', frameon = False)
plt.savefig(args.output_folder.rstrip("/") + '/distrib_compared' + prefix1 + '_vs_' + prefix2 + '.png')
print(f"Generated IBD segments length comparison plot")
    
