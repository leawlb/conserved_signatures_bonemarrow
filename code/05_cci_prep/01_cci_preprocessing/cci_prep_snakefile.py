#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/01_cci_preparation"

COLORS = config["base"] + config["metadata_paths"]["colors"]

LRDB_INP =  config["base"] + config["metadata_paths"]["lrdb_inp"]
LRDB_OUT =  config["base"] + config["metadata_paths"]["lrdb_out"]
ASSIGNMENT = config["base"] + config["metadata_paths"]["assignment"]
print(ASSIGNMENT)

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

VALUES = config["values"]["05_cci_prep"]

# from which colData slot to take identites (e.g. cluster, subcluster, cell type) as specified in config
COLNAME_IDENTITY = VALUES["colname_identity"]

#-------------------------------------------------------------------------------

targets = []
for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "/01_sprp/sce_" + s + "_" + a + "-01"]
    targets = targets + [OUTPUT_DAT + "/02_asnm/sce_" + s + "_" + a + "-02"]
    targets = targets + [OUTPUT_DAT + "/04_subs/sce_" + s + "_" + a + "-04"]
    targets = targets + [OUTPUT_DAT + "/05_down/sce_" + s + "_" + a + "-05"]

targets = targets + [LRDB_OUT]

targets = targets + [OUTPUT_REP + "/sce_downsampling_report.html"]

#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------


""" 
PREP SCE

The SCE objects need to be prepared before CCI calculation.
Certain conditions must be met for CCI calculation:

- cluster or cell type labels:
  - should be in colData slot called "Identity"
  - should not be able to grep each other (e.g. "B cell" and "pre-B cell")
  - should not contain ungreppable characters (+, (, ), or &)
  - should be factored (preferably in a reasonable sequence)

- RowData
  - rownames should contain Gene Symbols
  - a rowData slot called "ENSMUS_ID" should countain Gene IDs 

- Emitter/Receiver Assignment:
  - should be in colData slot called "Assignment"
  - should contain "emitter" for designated emitter cells
  - should contain "receiver" for designated receiver cells  
  
- Quality
  - all SCE objects to be compared should have comparable counts matrices
  - downsampled to lowest quality samples
"""

#-------------------------------------------------------------------------------

"""
Add Clusterlabels

- add labels from annotated SCE objects to age merges, store in "Identity" slot
- store ENSMUS IDs in "ENSMUS_ID" rowData slot

This step is individual for each SCE object, but after this step, 
following rules should be readily usable by adjusting assignment.txt as required
"""
rule add_clusterlabels:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/ages/sce_{species}_{age}-08",
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        sce_output = OUTPUT_DAT + "/01_sprp/sce_{species}_{age}-01"
    params:
        colname_identity = COLNAME_IDENTITY
    script:
        "scripts/01_add_clusterlabels.R" 

#-------------------------------------------------------------------------------

"""
Add Assignments

- add assignments for each Identity (emitter, receiver or remove) 
- from assignments.txt, adjust it manually as required
- unwanted cell types labeled "remove" are removed
- cell types with fewer than X cells are also removed, see config

- Identity is factorized according to Identity sequence in assignment.txt

!Because cell types with few cells are removed, not all CCI objects might later
contain the exact same cell types!
If all cell types should be contained in CCI objects, set min_cells to 0
"""

rule celltype_assignment:
    input:
        sce_input = rules.add_clusterlabels.output,
        assignment = ASSIGNMENT
    params:
        min_cells = VALUES["min_cells"]
    output:
        sce_output = OUTPUT_DAT + "/02_asnm/sce_{species}_{age}-02"
    script:
        "scripts/02_add_assignment.R" 
        
#-------------------------------------------------------------------------------

""" 
Subsetting

Not strictly required but useful for quicker computation and better downsampling

- SCE and LRDB objects are subset to intersecting genes 
- CellTalkDB "mouse_lr_pair.rds"" downloaded from 
https://github.com/ZJUFanLab/CellTalkDB/blob/master/database/mouse_lr_pair.rds
on Jul 6, 2023 at 2:50PM

"""

rule prepare_lrdb:
    input:
        sce_input_path = expand(rules.celltype_assignment.output, species = species, age = age),
        lrdb_input = LRDB_INP
    output:
        lrdb_output = LRDB_OUT # saved in metadata, not data
    script:
        "scripts/03_prepare_lrdb.R" 

rule subset_sce:
    input:
        sce_input = rules.celltype_assignment.output,
        lrdb_input = rules.prepare_lrdb.output
    output:
        sce_output = OUTPUT_DAT + "/04_subs/sce_{species}_{age}-04",
    script:
        "scripts/04_subset_sce.R" 
        
#-------------------------------------------------------------------------------

"""
Downsampling

- downsample all SCE objects to lowest quality sample to ensure comparability 
- downsampled counts are stored in a "downsampled" assay slot
- after subsetting so that unused genes don't cause quality differences
"""    

rule downsample:
    input:
        sce_input_path = expand(rules.subset_sce.output, species = species, age = age),
        assignment = ASSIGNMENT
    output:
        sce_output_path = expand(OUTPUT_DAT + "/05_down/sce_{species}_{age}-05", species = species, age = age)
    script:
        "scripts/05_downsample_ass.R" 

# make a report for downsampling QC
rule make_downsampling_report:
    input:
        sce_input_path =  expand(rules.downsample.output, species = species, age = age),
        lrdb_input = rules.prepare_lrdb.output
    output:
        OUTPUT_REP + "/sce_downsampling_report.html"
    params:
        colors_path = COLORS,
        sce_functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "sce_downsampling_report.Rmd"
        
