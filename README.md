# IRF4 Multiple Myeloma TRN
Code relevant to 2024 MM TRN IRF4/p300 manuscript

The repo is organized as following:

-ChIP processing (primarily in python), includes notebooks related to processing ChIP data files with relevant tools (including EnhancerPromoter, ROSE2, etc.). Also contained are log outputs from various tools used. This folder contains scripts and data folders for the public data processed, as well as the data processed specifically for this paper. The data folders includes the various intermediate steps for the ChIP data used. 

-ChIP data contains frequently referenced data within other notebooks. 

-RNA data (R), includes RNA data, notebooks used for LogFC calculation between WT and tested conditions, as well as some accompanying QC figures.

-Public data includes data that was collected from public resources

-figure n contains code used to generated figures observed within the manuscript. Generally base figures were made and then aesthetics were updated within Adobe Illustrator for figure generation. Figures generated in R are in the correct directory structure. Bam plots generated in python, are put within their relevant figure folder, but paths should be directed torwards the original ChIP processing folder. 

Some data is not included here (for example DepMap 22Q4 Chronos scores, or some of the MACS1.4 outputs) due to sizing being over the available sizing limits for github. 

12/2024 edit - following reviews, additional APMS was added to the manuscript. I did not have time while at Kronos to update these to each figure n, but APMS related figures can be found in "AP_MS" folder
