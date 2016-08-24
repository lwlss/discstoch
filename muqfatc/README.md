###### Python package \code{muqfatc} version 0.1

Link: https://github.com/lwlss/discstoch
Contact: h.herrmann2@gmail.com

#### General Usage 

This package provide an image analysis tool for capturing QFA (Addinall et al., 2011) cultures on a high-throughput level (i.e. muQA). Single lineage and pin population growth is catpured.

The package consist of four modules:

1. \code{imageanalysis} provides generic functions required for image analysis of yeast cells grown on solid agar. The two main functions are \code{makeBackground} for generating a smooth background of medium pixel intensity from the input folders containing microscopic image observations and \code{makeBorder} which draws a border around microscopic observations for better detection of cells growing at the side or off the side of the agar plate. 

2. \code{pingrowth} provides all of the functions required for capturing the total area of yeast on agar over time for each pin. Pin images are provided as folder inputs where each folder consists of a series of tiff images. 

3. \code{lineagegrowth} provides all the functions required for capturing single lineage time courses whereby each single cell present at the first time point is tracked over time. Merged colonies are not tracked. 

4. \code{resratemap} is an adhoc feature which after growth rate parameters of single lineages have been calculated from the obtained features can be specified as input. Colonies are colour-coded according to the computed rate and residuals. This feature was designed in order to be able to check whether image quality or growth rate patterns can be detected across the agar plate. 

#### Example Code 

All code I wrote which makes use of this specific package is provied in the Analyses/ImageAnalysis/ folder on GitHub along with a sample input folder R09C08 which contains the corresponding tiff images of the microscopic observations for that pin. 

 