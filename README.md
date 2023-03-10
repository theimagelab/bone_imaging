# Longitudinal intravital imaging code supplement.

This code accompanies the paper "Minimally invasive longitudinal intravital imaging of cellular dynamics in intact long bone"

Code adapted from [1]

### Instructions
1. Download the above R script (SurfDimReduc.R) and data directory (dimredclean_statistics) into the same folder on your system.
2. Open the R script in RStudio and set the working directory to the same directory as the R script. 
3. Install Seurat and dplyr R packages if you do not have them installed already.
4. Source the script.
Typical run time on a 2022 M1 Macbook Pro: <10 seconds.

### Expected output
1. Heatmaps as in Figure 5A 
2. Dimension Reduction plots as in Figure 5B.

### Dependencies
1. dplyr
2. Seurat

### Environment Information
This code has been tested using the following:
1. dplyr version 1.0.9
2. RStudio version 2022.07.2 
3. macOS Monterey v12.0
4. 2022 Apple M1 Macbook Pro

### References
1. Crainiciuc G, Palomino-Segura M, Molina-Moreno M, Sicilia J, Aragones DG, Li JLY, Madurga R, Adrover JM, Aroca-Crevillén A, Martin-Salamanca S, et al: Behavioural immune landscapes of inflammation. Nature 2022, 601:415-421.
