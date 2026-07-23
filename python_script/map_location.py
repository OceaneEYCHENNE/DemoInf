import argparse  
import tskit 
import pyslim
import random 
import numpy as np
import pandas as pd 
import matplotlib 
from matplotlib import pyplot as plt
from pathlib import Path

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to display the sampled individuals among all the individuals on the map")

# Add required arguments 
parser.add_argument("--input_file", "-i", type=str, required=True, help="Path for IBD file with coordinate ")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output folder")
parser.add_argument("--coord_sim_file", "-coord_sim", type=str, required=True, help="Path for simulated location file")
parser.add_argument("--coord_file", "-coord", type=str ,required=True, help="Path for sampled location file")
parser.add_argument("--map_file", "-map", type=str, required=True, help="Path for map file")
parser.add_argument("--time_MRCA", "-t", action="store_true", help="Map as a function of time since their MRCA")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD_coord_")[1]


### Load files ### 
locs_sim = pd.read_csv(args.coord_sim_file) 
locs = pd.read_csv(args.coord_file)

# Load a background map 
im = plt.imread(args.map_file)
fig, ax = plt.subplots()
ax.imshow(im)  															  


### Print points cloud plot ###
ax.scatter(locs_sim["X"], locs_sim["Y"], s=20, facecolors='none', edgecolors='grey', label='Simulated')    # simulated individuals in grey 
ax.scatter(locs["x"], locs["y"], s=20, c="red", label='Sampled')  									     # sampled individuals in red  
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_xlim(0, 100)                                                                       # Places the x-axis in ascending order (from 0 to 100)
ax.set_ylim(0, 100)																		  # Places the y-axis in ascending order (from 0 to 100)
ax.set_title('Location of sampled individuals among all simulated individuals')
legend = ax.legend(loc='upper right', frameon=True, title='Individuals', handletextpad=0.1, bbox_to_anchor=(1.02, 1.02))
legend.get_frame().set_facecolor('white')  
legend.get_frame().set_edgecolor('black')  


# Export map in a .png file 
plt.savefig(args.output_folder.rstrip("/") + '/Map_individuals_' + prefix + '.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"Generated map of individual locations")

### Map pairs of individuals sharing IBD segments as a function of time since their MRCA
if args.time_MRCA:	
	# Load IBD files 
	dfIBD = pd.read_csv(args.input_file)
	
	# Create fig 
	fig, ax = plt.subplots()
	ax.imshow(im) 
	
	# Tracing links between individuals sharing IBD segments of less than 10 generations
	for g in range(1, 11):  
		subset = dfIBD[dfIBD["Time"] == g]  											  # Filter individuals for g generation
		for _, row in subset.iterrows():  												  # Browse each line of the subset
			ax.plot([row['X1'], row['X2']], [row['Y1'], row['Y2']], color=plt.cm.autumn_r((10-g) / 10), linewidth=0.6)
	ax.scatter(locs["x"], locs["y"], s=20, c="red", zorder=3)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_title('IBD Connections Between Individuals')
	ax.set_xlim(0, 100)                                                                    
	ax.set_ylim(0, 100)	
	colors = [plt.cm.autumn_r((10-g) / 10) for g in range(1, 11)] 						  # Create a list of colors for each generations
	handles = [plt.Line2D([0], [0], color=color, linewidth=2) for color in colors]
	labels = [str(g) for g in range(1, 11)]
	ax.legend(handles=handles, labels=labels, title="Gen.", loc="lower right", frameon=True, title_fontsize=10, handletextpad=0.3)
	legend.get_frame().set_facecolor('white')  
	legend.get_frame().set_edgecolor('black')
	plt.savefig(args.output_folder.rstrip("/") + '/Map_individuals_time_MRCA_' + prefix + '.png', dpi=300, bbox_inches='tight')
	plt.close()
	print(f"Generated map of individual locations with time since their MRCA")

        
        
