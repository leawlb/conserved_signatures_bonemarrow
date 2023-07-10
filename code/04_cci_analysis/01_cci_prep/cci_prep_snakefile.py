#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_prep/"
OUTPUT_REP = OUTPUT_BASE + "/reports/04_cci_prep/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

targets = []

for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "01_sprp/sce_" + s + "_" + a + "-01"]
    targets = targets + [OUTPUT_DAT + "02_asnm/sce_" + s + "_" + a + "-02"]
    targets = targets + [OUTPUT_DAT + "04_subs/sce_" + s + "_" + a + "-04"]
    targets = targets + [OUTPUT_DAT + "05_down/sce_" + s + "_" + a + "-05"]
    targets = targets + [config["metadata"]["path"] + "/CCI/lrdb"]

targets = targets + [OUTPUT_REP + "sce_downsampling_report.html"]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
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
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/02_clst/louvn_clust/sce_hsc-test",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/02_clst/louvn_clust/sce_str-test"
    output:
        sce_output = OUTPUT_DAT + "01_sprp/sce_{species}_{age}-01"
    script:
        "scripts/01_add_clusterlabels.R" 


"""
Add Assignments

- add assignments for each Identity (emitter, receiver or remove) 
- from assignments.txt, adjust it manually as required
- unwanted cell types labeled "remove" are removed
- cell types with fewer than X cells are also removed, see config

-Identity is factorized according to Identity sequence in assignment.txt

!Because cell types with few cells are removed, all CCI objects might later
not contain the exact same cell types!
If all cell types should be contained in CCI objects, set min_cells to 0
"""

rule celltype_assignment:
    input:
        sce_input = rules.add_clusterlabels.output,
        assignment = "assignment.txt"
    params:
        min_cells = config["values"]["cci_prep"]["min_cells"]
    output:
        sce_output = OUTPUT_DAT + "02_asnm/sce_{species}_{age}-02"
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
        lrdb_input = config["metadata"]["path"] + "/CCI/mouse_lr_pair.rds"
    output:
        lrdb_output = config["metadata"]["path"] + "/CCI/lrdb"
    script:
        "scripts/03_prepare_lrdb.R" 

rule subset_sce:
    input:
        sce_input = rules.celltype_assignment.output,
        lrdb_input = rules.prepare_lrdb.output
    output:
        sce_output = OUTPUT_DAT + "04_subs/sce_{species}_{age}-04",
    script:
        "scripts/04_subset_sce.R" 
        
#-------------------------------------------------------------------------------

"""
Downsample

- downsample all SCE objects to lowest quality sample to ensure comparability 
- downsampled counts are stored in a "downsampled" assay slot
- after subsetting so that unused genes don't cause quality differences
"""    

rule downsample:
    input:
        sce_input_path = expand(rules.subset_sce.output, species = species, age = age)
    output:
        sce_output_path = expand(OUTPUT_DAT + "05_down/sce_{species}_{age}-05", species = species, age = age)
    script:
        "scripts/05_downsample.R" 

# make a report for downsampling QC
rule make_downsampling_reports:
    input:
        sce_input_path =  expand(rules.downsample.output, species = species, age = age),
        lrdb_input = rules.prepare_lrdb.output
    output:
        OUTPUT_REP + "sce_downsampling_report.html"
    script:
        "sce_downsampling_report.Rmd"
        
