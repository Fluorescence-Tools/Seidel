
# *Seidel* software for FRET-nanoscopy analysis  
## functionality  
1) reading photon streams from picoquant .ptu files and transforming them into xy-lifetime images
2) Bayesian framework to optimally determine the number and location of emitters in an image snapshot
3) Forming individual emitters into FRET pairs and calculating FRET parameters   
4) Bayesian framework for fitting 1D histograms with non-centered Chi and Gaussian distributions.
5) Molecular assembly particle averaging based on coarse alignment and Prokrustes analysis

## System requirements
Software consists of python code with a c++ backend. The latter is responsible for computation intensive tasks of reading ptu files, creating images and performing gaussian fits. Python software was build on python 3.7 using libraries available in the *anaconda* standard cluster of libraries. c++ code was compiled into a .pyd file that requires windows 10 standard dll's to run.
* python 3.7
* anaconda libraries
* windows 10  
Software was tested by several users under these settings, no non-standard hardware is required.

## installation guide
1) clone the /FRC/Code into a local directory
2) include the path to your local copy of /FRC/Code in your PYTHONPATH  
3) Seidel depends on the *lmfit scikit-image numpy scipy matplotlib tiffile* packages included in the metapackage anaconda. The package *rmsd* must be downloaded separately using
    pip install rmsd
In case of any problems, please contact voort@hhu.de. Installation time should be less than a minute.

## Demo notebook & example data
For a step-by-step instruction on how to use *Seidel* on FRET nanoscopy data, refer to [this](https://github.com/Fluorescence-Tools/FRC/blob/master/templates/20211103_template_v3.ipynb) notebook. The notebook shows expected results for comparison. Run-time is ~1 second per ptu file.
A small dataset of ptu files is provided. As well as a slightly larger processed dataset that can be used for further analysis. To obtain a complete FRET-nanoscopy dataset of several Gigabytes, please contact voort@hhu.de. Similarly, to reproduce published results obtained using *Seidel*, please request the relevent dataset and accompanying analysis notebook.

## dedication
The *Seidel* software is dedicated to Prof. C.A.M. Seidel from the Institute of Molecular Physcial Chemistry at the HHU university, the PhD supervisor of the author.
