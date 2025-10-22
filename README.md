# DiverseCollections
A repository to house the scripts, appropriate data, and other associated content for the paper titled: Diverse collection of wild Drosophila melanogaster reveals the evolutionary history of a persistent viral infection


Here is an overview of how scripts and data files in this repository map to figures and tables in the paper:


### Figure 1

Maps were created using [create_sample_map.R](analyses/scripts/create_sample_map.R) using data files referenced therein.  In general, scripts refer to file paths relative to the root of this repository.  This can be created by setting the working directory to Project Directory in RStudio.

qPCR data panels were created using [qPCR_data.R](analyses/scripts/qPCR_data.R) using input datafiles referenced therein.  The map images and qPCR data panel;s were combined using [Affinity Designer](figures/affinity_designer/Figure_1.afdesign) and output as [a PDF](figures/Figure_1.pdf).

Supplemental Tables 2 and 3 were also output from the qPCR_data.R script




