# "The landscape of actionable genomic alterations in cell-free circulating tumor DNA from 21,807 advanced cancer patients" published in Clinical Cancer Research, DOI: 10.1158/1078-0432.CCR-17-3837

This repository contains scripts and data for generating the figure panels from the Zill et al 21,807 patient ctDNA landscape study.

Prior to running all analyses involving comparison to TCGA data, the TCGA variants were downloaded from www.cbioportal.org using the get_cbio_data.r Rscript.  

Note: that the TCGA data for the study were downloaded on 4/25/16, using the get_cbio_data.r script. Doing this now will produce a non-equivalent analysis, as the TCGA data have been updated in cBio portal since 4/25/16.

##############################
# regenerating figure panels #
##############################

To regenerate figure panels, do the following from the command line:

1) Clone the repository to you local instance

	> git clone https://github.com/guardant/Zill_2018.git

2) change directories to Zill_2018# all additional commands (Rscript calls) must be executed from within the main package directory
3) Call "Rscript R/ctdna_landscape_analysis.R"
4) All figure panels and tables will be written to /output

#######################################################
# Relationship between scripts and manuscript figures #
#######################################################

Helper scripts:
---------------
get_cbio_data.r		#get TCGA data from cbio portal.
scatterBarOZ.R 		#used by pancan_clonal_frxn.r

Figures code to generate individual figures is commented in ctdna_landscape_analysis.R

##############
# Input data #
##############
	
Data for generating the figures is located in data_tables_for_figures directory.
