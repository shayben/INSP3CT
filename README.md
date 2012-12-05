INSP3CT
=======

Interpolation and Statistical Proximity of 3C Tables

Version 1.0
Copyright 2012 Shay Ben Elazar ©
--------------------------------------------------------------------------------
Provided as supplementary alongside - 
"Spatial localization of co-regulated genes exceeds genomic gene clustering in 

the S. cerevisiae genome"
By Shay Ben-Elazar, Zohar Yakhini, and Itai Yanai. Accepted for publication to 

NAR (Dec 2012).
--------------------------------------------------------------------------------

Description:
-------------
This program can be used to take a sparse contact dataset and interpolate it to 

a sampled matrix at an arbitrary resolution.
On top of the resulting matrix, any set of annotation of loci is then inspected 

for co-localization both in 3d and 1d.
The result can then optionally be displayed in a figure similar to Figure 3 in 

the above paper.

Requirements:
--------------
It is recommended to run this on a machine with at least 8 GB of RAM, depending 

on the size of the dataset.
Multiple cores can and should be utilized to speed up analysis.
Our setup (for this particular dataset) was a machine with 192 GB RAM and 12 

dedicated cores.
Basic Matlab knowhow is a prerequisite.

Input and parameters:
----------------------
Edit the main pipeline file, INSP3CT.m, and make sure to set the following 

parameters.
Required input files:
* PARAMETERS.interfrequencies - filename of the inter-chromosomal frequencies 

tab-delimited in the format with the following header:   chr1	locus1	chr2	

locus2	frequency
* PARAMETERS.intrafrequencies - filename of the intra-chromosomal frequencies 

tab-delimited in the same format, where chr1=chr2.
* PARAMETERS.ChromosomeSizes - An array holding the number of bases in each 

chromosome in the genome.
	
An example of a sufficient file-set can be found with the supplementary 

materials (YeastDatasetExample.zip). The corresponding parameters are commented 

out in the code.

Additional parameters:
* PARAMETERS.cores - number of cores in the machine used for processing. More 

cores is faster but requires more RAM as well.
* PARAMETERS.binningresolution - the number of nucleotides which are considered 

as a locus. Smaller value -> larger dataset and slower to compute but provides 

more accuracy in the results.
* PARAMETERS.analysis_shuffle - set to true for a random control.
* PARAMETERS.analysis_cyclic - set to true for a cyclic permutation control.

You also need to manually create the following variables before running any co-

localization analysis:
* Loci - A vector of k bin indices containing points of interest.
* AnnotationData - A matrix of size kXn where each row is a point of interest 

and each column an annotation. The (i,j)-th entry is the number of elements 

annotated as j which reside in bin Loci[i].
There is a commented-out example which randomly generates loci and annotations 

used for testing.

Usage:
-------
Main file is INSP3CT.m.
* The first code block, labeled "INSP3CT pipeline", loads the parameters and 

asserts their validity.
* The second code block, labeled "Load and parse datasets", Prepares the data 

structures used to generate the interpolated contact frequency matrix.
* The third code block, labeled "Generate interpolated and binned frequencies", 

runs the interpolation step.
The next two blocks can be used if Loci and AnnotationData are set:
* The forth code block, labeled "Analysis on annotation" given an annotation on 

the binned genome, analyzes 3D vs 1D co-localizations using mHG (see Eden, E., 

Lipson, D., Yogev, S. and Yakhini, Z. (2007) Discovering motifs in ranked lists 

of DNA sequences. PLoS Comput Biol, 3, e39 and Eden, E., Navon, R., Steinfeld, 

I., Lipson, D. and Yakhini, Z. (2009) GOrilla: a tool for discovery and 

visualization of enriched GO terms in ranked gene lists. BMC Bioinformatics, 10, 

48.).
* The fifth code block, labeled "Show a result figure for each annotation", 

plots each annotation showing the resulting landscape of enrichment. Pressing 

any keyboard key will continue to show the next annotation.

Output:
--------
obj.Interpolated - cell array containing interpolated chromosomal contact 

frequencies binned to the requested resolution.
If Loci and AnnotationData are set:
The resulting enrichment comparisons are a set of figures for the respective 

annotations.
corrected3d - Corrected mHG (partition limited) score for 3D enrichment around 

loci for each annotation.
corrected1d - Corrected mHG (partition limited) score for 1D enrichment around 

loci for each annotation.
allbestspheresizes - mHG threshold on 3D ordered loci.
allbestintsizes - mHG threshold on 1D ordered loci.

Disclaimer:
------------
These scripts are provided as is and under the GPL license.
Much of the code is *NOT* optimized since it has undergone many iterations and 

was not specially tailored for the purposes in this code. Please make sure to 

read the instructions and paper thoroughly before contacting the authors.
Importantly, please note that the analysis in the paper was performed on a 

haploid genome and there are repercussions for switching to a diploid genome in 

the interpretation of the input data which were not addressed in this code or 

paper.

Copyright:
----------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact:
---------
For questions regarding this program e-mail shayben@cs.technion.ac.il
Further inqueries can be directed to Itai Yanai, web: 

http://yanailab.technion.ac.il/ email: yanai@technion.ac.il
