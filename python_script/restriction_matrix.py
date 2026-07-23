import pandas as pd
import glob
import os 

### Load SLiM matrix ### 
slim_matrix = pd.read_csv("/mnt/DATA1/MAP/map_100x100_grad_param_value.csv")

# Spatial subsample -> remove grid boundaries
slim_sub = slim_matrix[
    slim_matrix['x_km'].between(1, 98) & 
    slim_matrix['y_km'].between(1, 98)
]

### Get all MAPS files 
maps_files = glob.glob("/mnt/DATA1/MAPS_RESULTS_KM/grad_N_m1.0_N5.0_theo_K200_t201/all_MAPS_Ne*.csv")

### Filter + Merge with SLiM values 
results = []

for filepath in maps_files:

    maps_matrix = pd.read_csv(filepath)

    basename = os.path.basename(filepath)
    prefix = basename.replace("MAPS_Ne_", "").replace(".csv", "")

    maps_sub = maps_matrix[
        maps_matrix['x_km'].between(1, 98) &
        maps_matrix['y_km'].between(1, 98)
    ]

    # Merge MAPS values with thoritical values of SLiM
    merged = pd.merge(
        slim_sub[['x_km', 'y_km', 'true_value']],
        maps_sub[['x_km', 'y_km', 'x', 'y', 'value']],
        on=['x_km', 'y_km'],
        how='inner'
    )

    outpath = f"/mnt/DATA1/MAPS_RESULTS_KM/grad_N_m1.0_N5.0_theo_K200_t201/MAPS_Ne_{prefix}.csv"
    
    merged.to_csv(outpath, index=False)



