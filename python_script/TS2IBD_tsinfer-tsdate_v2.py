import argparse  
import tskit 
import msprime
import random 
import numpy as np
import pandas as pd 
import re
from pathlib import Path
import json


### Define arguments ### 
# Initialize parser
parser = argparse.ArgumentParser(description="Script to obtain IBD segments file from a tree-sequence file generated from tsinfer/tsdate")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str,required=True, help="Path for input tree-sequence")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output ")
parser.add_argument("--haplotypes", "-hap",type=int, default=200, help="Number of haplotypes sampled : 2 times number of individuals sampled")
parser.add_argument("--max_time", "-t_max", type=int, default=None, help="Maximal ages for IBD segments recorded")
parser.add_argument("--min_span", "-l_min",type=int, default=None, help="Minimal length for IBD segments recorded")
parser.add_argument("--location", "-loc",  type=str, default=None, help="CSV file with node coordinates")

args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
prefix = input_path.name.split(".trees")[0]

# Modify the prefix for integrate the information of haplotype and max time 
if "_rep" in prefix:
    before_rep, after_rep = prefix.split("_rep", 1)
    prefix = f"{before_rep}_K{args.haplotypes}_t{args.max_time}_rep{after_rep}"
else:
    # If no "_rep", simply add at the end
    prefix = f"{prefix}_K{haplotypes}_t{max_time}"

### Import tree-sequence of an SLiM simulation ###
ts = tskit.load(args.input_file)
print(f"Loaded a simulated tree sequence of {ts.num_samples} haploid genomes")
print(f"containing {ts.num_sites} variable sites over {ts.sequence_length} base pairs")

### Obtain IBD segments ### 
# Stock IBD segments
segments = ts.ibd_segments(store_pairs=True, store_segments=True, max_time = args.max_time, min_span = args.min_span)

# Define IBD segment stitching function
def stitch_ibd(ts):
	# Create numpy array empty for stock values 
	Pair = []																      		  # create a list of haplotypes pairs 
	Left = np.array([])																      # create a numpy array of start positions of IBD segments  
	Right = np.array([])																  # create a numpy array of end positions of IBD segments
	Node = np.array([])																      # create a numpy array of ID node of IBD segments
	Time = np.array([])																      # create a numpy array of age of the nodes
	
	# Traverses each element by each haplotypes pairs (k) and its segments (v)
	for k, v in segments.items():
		kLeft = v.left 																      # store the start positions of segments 
		kRight = v.right 																  # store the end positions of segments
		kNode = v.node															  		  # store ID node of IBD segments
   		
   		# Order segments by ascending order of start position
		order = kLeft.argsort()												  		      # return the sort indices  
		l = kLeft[order]														  		  # order the start position ordered
		n = kNode[order]														  		  #	order the node ID ordered
		r = kRight[order]														  		  # order the end positions ordered
    	
    	# Define convolution function : detect breaks in the IDs of ancestral nodes in order to group together consecutive segments belonging to the same ancestor and cut out where the ancestor changes
		node = (np.convolve(n + (n[0] == 0), [-1, 1]) != 0)[:-1]					      # convolution function that detects discontinuities in ID nodes of two adjacents iIBD segments
		discontinuities = np.concatenate(([True], l[1:] != r[:-1]))						  # detect discontinuities in left and right positions of two adjacents IBD segements 
		i = np.where(node | discontinuities)[0]											  # find indices of dicontinuities in ID nodes and/or positions start/end 
        
    	# Associate each element with the pair k 
		kPair = [k, ]*len(i)														      # create a list of the pair k repeted 
		kLeft = l[i]														  		      # select start positions at i indices 
		kRight = np.append(r[i[1:]-1], r[len(r)-1])										  # adjust end positions 
		kNode = n[i]   														  		 	  # select ID nodes at i indices 
    	
    	# Add values at existing numpy array 
		Pair = Pair + kPair														  		  # add the new pair at the list 
		Left = np.append(Left, kLeft)													  # add start positions 
		Right = np.append(Right, kRight)											      # add end positions
		Node = np.append(Node, kNode)													  # add ID nodes 
		Time = np.append(Time, ts.tables.nodes.time[kNode])								  # add age of the nodes 
	
	# Calculate the length of IBD segments 
	Length = np.subtract(Right,Left)	
	
	# Clean the name of Pair 
	Pair = [(int(p[0]), int(p[1])) for p in Pair]
    
	# Create a dataframe containing all the informations of the IBD segments 
	dfIBD = pd.DataFrame(columns = ['Pair', 'Hap1', 'Hap2', 'Ind1', 'Ind2', 'Left', 'Right', 'Length', 'Node', 'Time'])  
	dfIBD['Pair'] = Pair														  		  # fill column Pair with haplotype pairs                                  														
	dfIBD['Left'] = Left														  		  # fill column Left with start positions 
	dfIBD['Right'] = Right														  		  # fill column Right with end positions 
	dfIBD['Length'] = Length														  	  # fill column Length with IBD length  
	dfIBD['Node'] = Node													  	  	      # fill column Node with Age of IBD 
	dfIBD['Time'] = Time														  		  # fill column Time with IBD age 
	
	return dfIBD
	
# Stitch IBD segments 
dfIBD = stitch_ibd(ts)

# Create new columns containing additional information on haplotypes pairs  
dfIBD[['Hap1', 'Hap2']] = pd.DataFrame(dfIBD['Pair'].tolist(), index=dfIBD.index)    # create columns containing haplotypes IDs of each pair 

# Create new columns containging additional information on individuals pairs 
node_to_ind = {node_id: ts.node(node_id).individual for node_id in ts.samples()}          # build a map from nodes to individuals
dfIBD["Ind1"] = dfIBD["Hap1"].map(node_to_ind)                                            # apply this map for record the individual ID of the halotype1
dfIBD["Ind2"] = dfIBD["Hap2"].map(node_to_ind)                                            # apply this map for record the individual ID of the halotype2

# Keep only IBD segments of individuals pair which are alive at time 0 
focal_ind_ids = set(ind.id for ind in ts.individuals() if ind.time == 0)
dfIBD = dfIBD[dfIBD['Ind1'].isin(focal_ind_ids) & dfIBD['Ind2'].isin(focal_ind_ids)]

# Return the number of IBD segments recorded 
num_ibd = len(dfIBD.Pair) 
print(f"Recorded {num_ibd} IBD segments")

# Export dataframe into .csv in output folder 
dfIBD.to_csv(args.output_folder.rstrip("/") + '/IBD_' + prefix + '.csv', index=False)


### Add individuals coordinates ### 
if args.location:	
	# Import coordinates file 
	locs_all_csv = pd.read_csv(args.location, sep = ",")
	locs_all_csv.to_csv(args.output_folder.rstrip("/") + '/Coord_all_' + prefix + '.csv', index=False)
	
	# Recover indivudals coordinates of sampled individuals
	focal_inds = [ind for ind in ts.individuals()]
	ids = [ind.id for ind in focal_inds]
	vcf_ids = []
	for ind in focal_inds:
	    try:
	        meta = json.loads(ind.metadata.decode("utf-8"))
	        vcf_ids.append(meta.get("uid", str(ind.id)))
	    except:
	        vcf_ids.append(str(ind.id))
	        
	locs_all_csv["id"] = locs_all_csv["id"].astype(str)
	locs_csv = locs_all_csv[locs_all_csv["id"].isin(vcf_ids)].set_index("id").loc[vcf_ids]

	# Save coordinates into a dataframe 												
	dfLOC = pd.DataFrame(columns=['Ind_ts', 'Ind_vcf', 'X', 'Y'])
	dfLOC['Ind_ts'] = ids  									  			      # IDs of alive individuals at generation 0 
	dfLOC['Ind_vcf'] = vcf_ids 
	dfLOC['X'] = locs_csv["x"].values  															  # Coordinates X
	dfLOC['Y'] = locs_csv["y"].values 															      # Coordinates Y
	dfLOC = dfLOC.set_index('Ind_ts')

	Num_Ind = len(dfLOC.Ind_vcf)  															  # Total number of recorded individuals
	print(f"Recorded {Num_Ind} individual locations")
	
	### Add column coordinates in the dataframe ###
	# Extract the pairs of haplotype identifiers 
	hap_pair = np.array(list(zip(dfIBD["Hap1"], dfIBD["Hap2"])))

	# Get the corresponding pair of individuals
	ind_pair = np.column_stack((dfIBD["Ind1"].to_numpy(), dfIBD["Ind2"].to_numpy()))

	# Get the coordinates of the first individual in pair
	X1 = dfLOC.loc[ind_pair[:, 0]]['X'].values
	Y1 = dfLOC.loc[ind_pair[:, 0]]['Y'].values
    
	# Get the coordinates of the second individual in pair
	X2 = dfLOC.loc[ind_pair[:, 1]]['X'].values
	Y2 = dfLOC.loc[ind_pair[:, 1]]['Y'].values

	# Add coordinates to the IBD dataframe
	dfIBD[['X1', 'X2', 'Y1', 'Y2']] = np.column_stack((X1, X2, Y1, Y2))

	# Export dataframe into a .csv file in output folder 
	dfIBD.to_csv(args.output_folder.rstrip("/") + '/IBD_coord_' + prefix + '.csv', index=False)
	dfLOC.to_csv(args.output_folder.rstrip("/") + '/Ind_Coord_' + prefix + '.csv', index=True)
