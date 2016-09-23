## Python package `muqfatc` version 0.1

Link: https://github.com/lwlss/discstoch ; https://pypi.python.org/pypi/muqfatc/0.1

Contact: h.herrmann2@gmail.com

This package is now on PyPI and can be downloaded using `pip install muqfatc`.

#### General Usage 

This package provides an image analysis tool for capturing QFA ([Addinall et al., 2011](http://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1001362)) cultures on a high-throughput level we call muQFA. Single lineage and pin population growth is captured.

The package consist of four modules:

1. `imageanalysis` provides generic functions required for image analysis of yeast cells grown on solid agar. For example, this modules includes functions for generating a smooth background of medium pixel intensity from the input folders containing microscopic image observations and for drawing a border around microscopic observations for better detection of cells growing at the side or off the side of the agar plate.

2. `pingrowth` provides all of the functions required for capturing the total area of yeast on agar over time for each pin. Pin images should be provided as folder inputs where each folder consists of a series of tiff images. 

3. `lineagegrowth` provides all the functions required for capturing single lineage time courses whereby each single cell present at the first time point is tracked over time. Merged colonies are not tracked. Again, pin images should be provided as folder inputs where each folder consists of a series of tiff images. Time course images for all colonies along with area estimate growth curves are returned.

4. `resratemap` is an ad-hoc feature whereby inferred growth rate parameters of single lineages and the corresponding residuals of model fitting can be specified as input. Colonies colour-coded according to the computed rate and residuals are returned. This feature was designed in order to be able to check whether image quality and growth rate patterns can be detected within a pin and across the agar plate. 

#### Example Code 

All code I wrote for my MSc Dissertation analyses which makes use of this specific package is provied in the Analyses/ImageAnalysis/ folder on GitHub along with a sample input folder R09C08 which contains the corresponding tiff images of the microscopic observations for that pin. This example data has not yet been integrated into the `muqfatc` package. 

#### Dependencies

This package relies on OpenCV (version 3.0.0) which is currently not installed as a dependency when downloading this package. OpenCV (cv2; version 3.0.0) must be manually installed when using this package! 

 
