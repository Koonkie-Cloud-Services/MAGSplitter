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
To use this script, set main.py as an executable file.  Following, pass the file locations for the input files as arguments.
A "results" folder will be created in the directory of MAGSplitter, and the output files will be stored there.
```
options:
  -h, --help            show this help message and exit
  -pf PF_FILE, --pf_file PF_FILE
                        metapathways ePGDB output file location (typically ptools/0.pf)
  -orf ORF_MAPPING, --orf_mapping ORF_MAPPING
                        metapathways duplicate ORF mapping file location (typically results/annotation_table/orf_map.txt)
  -contig ORF_CONTIG_MAP_FILE, --orf_contig_map_file ORF_CONTIG_MAP_FILE
                        metapathways orf to contig mapping file location (typically results/annotation_table/<samplename>.ORF_annotation_table.txt)
  -mag CONTIG_MAG_MAP, --contig_mag_map CONTIG_MAG_MAP
                        wgs pipeline contig contig to mag mapping file location (typically binning/results/greedy/config_info.tsv
```

```
example usage:
python main.py -pf example/0.pf -orf example/orf_map.txt -contig example/GAPP-5498e568-6918-4000-b27e-dbeff35eeee7.ORF_annotation_table.txt -mag example/contig_info.tsv
```

