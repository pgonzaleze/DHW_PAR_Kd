# Predicting cold-water bleaching in corals 

#### Repository with the R and python/jupyter notebook scripts, and data frame regarding the coral bleaching caused by cold sea surface temperature paper published on [MEPS doi.org/10.3354/meps13336](https://doi.org/10.3354/meps13336)

## Instructions

1.- Download environmental data in NetCDF format. The data is publicly available at:
  * [SST](https://coralreefwatch.noaa.gov/product/5km/index.php)
  * [PAR](https://oceancolor.gsfc.nasa.gov/l3/)
  * [Kd490](https://oceancolor.gsfc.nasa.gov/l3/) 

2.- Use any 'site_extraction_values*.ipynb' file to extract the desired environmental data for specified sites. After running the code, the output is saved in '.h5' format. You need the '*.h5' file to run the next step.

3.- Use the 'df_*-DCW.ipynb' file to compute the metrics. e.g DCW, ColdSpots, etc. The output is a dataframe/table

4.- Merge all the data-frames and save them as .csv 

5.- Use the 'logistic_regression.R' to perform the statistical tets *

*Notice the 'DCW_PAR_Kd.csv' file has the condensed data so that the 'logistic_regression.R' can be  used directly 

