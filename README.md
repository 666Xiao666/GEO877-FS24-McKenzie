# GEO877-FS24-McKenzie
# General information

Welcome to our GitHub!
Here you can find information about our code.

# 1-data_processing.ipynb
This file contains the preprocessing of our data. The first code block imports required modules, the second is responsible for data cleaning and the third creates the merged files for the four years together.
These codes to not have to be run, as all the data can already be found preprepared in the folder "data/cleaned_data/".

# 2-hexagon_analysis.ipynb
Here the code which generates our hexagon grid and the code which finds the intersections can be found. As well as calculation and map creation for results and analysis.

# 3-Voronoi_analysis.ipynb
Here the code which generates our Voronoi polygons can be found. Other than that the code is more or less identical to the code in 2-hexagon_analysis.ipynb, but small adjustments that everything works with Voronoi have been made.

# data
In the data folder the data we used is stored.

## original_data
This folder contains the original bike_share data downloaded from the website of cogo bike share.

## cleaned_data
This folder contains the data prepared with 1-data_processing.ipynb.

## Corporate_Boundary.geojson
This file contains a polygon with the boundaries of the city of Columbus.

# old
In this folder we archived files that were created but not used in the end or that were replaced by other code. This files are not needed for our algorithm anymore.