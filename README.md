# MAGSplitter
Script to split MAG files into individual MAG files.

### Background Information
***
Currently (August 2022), metapathways outputs are outputted in the form of ePGDBs readable for pathway tools
. In other words, the outputs signify entire metagenomes, rather than individual MAGs. If the end user were to 
want to look at metabolic pathway information for a single MAG, they would need to redo the metapathways pipeline
for each single MAG. This is inefficient in terms of time and resources. Therefore, this tool alleviates this 
inefficiency by taking the outputs from the WGS pipeline and the metapathways pipeline, and splitting the metabolic
ORF information into individual MAGs.

However, currently the file for ORF metabolic reactions within the ePGDB file has duplicate reactions previously removed.  THis is because 
the metapathways piepeline automatically removes duiplicate pathways since frequency information is not utilized in pathway tools.  Therefore, a step to bring these
ORFs are required beforehand.

Once the removed ORF reactions are added to our ePGDB file, the ORFs are then mapped to their respective metagenome contigs and MAGs.;
Following, the reactions are split into individual MAGs, and new 0.pf files are created for each MAG.


### Approach
***
* The ORF metabolic information outputted by metapathways is converted from PGDB format to a dataframe
* A list of missing ORFs are used to de-remove the removal of ORFs from the ePGDB file.
* A mapping file between ORFs and metagenome contigs are converted into a dataframe
* A mapping file between metagenome contigs and MAGs are converted into a dataframe
* Mapping files are then used to group the ORFs by MAGs.

### Inputs
***

* 0.pf: pathway tools initial input from metapathways outputs
(found in /ptools/)

* <sample_name>.ORF_annotation_tabel.txt: map of ORFS and contigs for an environment metagenome, outputted by metapathways 
(found in /results/annotation_table/)

* orf_map.txt: map of duplicated ORFS within a metagenome, with the first ORF being the intact ORF, 
and the rest being the duplicates
(found in /results/annotation_table/)

* config_info.tsv: map of contigs to MAGs created through the WGS pipeline binning process.  
(found in /binning/results/greedy)

### Usage 
***
MAGSplitter is meant to be used after the WGS and the metapathways pipeline.  
To use this script, open the "main.py" file with your favorite IDE.  Change the variables "pf_file",
"0rf_map_file", "orf_contig_map_file", and "contig_mag_map_file" into the appropriate paths for your file locations.  
Run "main.py", and the folders will be located in the /results folder.
