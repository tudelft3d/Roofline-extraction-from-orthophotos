# Unsupervised Roofline Extraction from True Orthophotos
This is the implementation part of the roofline extraction of the paper: [*Unsupervised Roofline Extraction from True Orthophotos*](https://arxiv.org/abs/2310.01067). 
First, we use a true orthophoto as input. 
Next, we perform line extraction to partition the building footprint, which generates separate roof parts. 
We then utilize a dense point cloud to extrude the partition results and create a LoD2 building model. For generating LoD2 building models, please refer to the repository [gfp-building-reconstruction](https://github.com/geoflow3d/gfp-building-reconstruction).

<div align="center">    
<img src="images/pipeline_all.jpg" width="800px" />
</div>

## Citation

If you use it in a scientific work, we kindly ask you to cite it:

<div class="filteredelement"><strong>Unsupervised Roofline Extraction from True Orthophotos for LoD2 Building Model Reconstruction</strong>. Weixiao Gao, Ravi Peters, and Jantien Stoter. <em> The 18th 3D GeoInfo 2023 </em>. <br/><a href="https://arxiv.org/abs/2310.01067"><i class="fas fa-external-link-alt"></i> PDF</a> <a href="#myref" data-toggle="collapse"><i class="fas fa-caret-square-down"></i> BibTeX</a> <div id="myref" class="collapse" tabindex="-1"><pre class="bibtex">@article{sum2021,
author = {Weixiao Gao, Ravi Peters, and Jantien Stoter},
title = {Unsupervised Roofline Extraction from True Orthophotos for LoD2 Building Model Reconstruction},
journal = {Lecture Notes in Geoinformation and Cartography (LNG&C) series},
year={2023},
publisher = {Springer},
}
</pre></div></div>


## Requirements 

*1.* Install all required Python packages as follows
```
pip install numpy shapely rtree tqdm joblib GDAL
```
The code was tested on numpy 1.24.2, shapely 2.0.1, rtree 1.0.1, tqdm 4.65.0, joblib 1.2.0, GDAL 3.6.2. 

*2.* Install GCC, Boost (1.63.0 or newer), Eigen3, CGAL, and OpenCV in Conda:
```
conda install -c conda-forge gcc=12.1.0; conda install -c anaconda boost; conda install -c omnia eigen3; conda install eigen; conda install -c conda-forge cgal; conda install -c conda-forge opencv
```
The code was tested on Eigen 3.3.7, CGAL 5.5.2, OpenCV 4.7.0.

*3.* Compile the ```libkinetic_partition.``` libraries:
```
CONDAENV=YOUR_CONDA_ENVIRONMENT_LOCATION
cd kinetic_partition
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DPYTHON_LIBRARY=$CONDAENV/lib/libpython3.10.so -DPYTHON_INCLUDE_DIR=$CONDAENV/include/python3.10 -DBOOST_INCLUDE_DIR=$CONDAENV/include -DEIGEN3_INCLUDE_DIR=$CONDAENV/include/eigen3 -DGMP_INCLUDE_DIR=$CONDAENV/include -DGMP_LIBRARY_RELEASE=$CONDAENV/lib/libgmp.so -DMPFR_INCLUDE_DIR=$CONDAENV/include -DMPFR_LIBRARIES=$CONDAENV/lib/libmpfr.so -DGMPXX_INCLUDE_DIR=$CONDAENV/include -DGMPXX_LIBRARIES=$CONDAENV/lib/libgmpxx.so -DCGAL_DIR=$CONDAENV/lib/cmake/CGAL -DOpenCV_DIR=$CONDAENV/lib/cmake/opencv4
make
```
The code was tested on Ubuntu 22.04 with Python 3.10.

### Data organization
The data in the folder is organized as follows:
```
data/img/*.tif    #orthophotos
data/poly/*.gpkg  #building footprints
...
```
The rest of the folders will be created automatically and the output rooflines will be stored in ```data/rooflines/*.gpkg ```.

## Running the code
```
source activate YOUR_CONDA_ENVIRONMENT
python3 main.py --ROOT_PATH=YOUR_DATA_LOCATION
```

## License
This implementation is free software; you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation; either version 3
of the License or (at your option) any later version. The full text of the license can be
found in the accompanying 'License' file.

If you have any questions, comments, or suggestions, please contact me at <i>gaoweixiaocuhk@gmail.com</i>

[<b><i>Weixiao GAO</i></b>](https://3d.bk.tudelft.nl/weixiao/)

May. 1st, 2023
