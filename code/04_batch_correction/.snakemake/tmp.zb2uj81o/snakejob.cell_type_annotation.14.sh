#!/bin/sh
# properties = {"type": "single", "rule": "cell_type_annotation", "local": false, "input": ["/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge/sce_mspr-06"], "output": ["/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_ctyp/sce_mspr-07"], "wildcards": {"species": "mspr"}, "params": {"ref_baccin_sce": "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_baccin_sce", "ref_dahlin_sce": "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_dahlin_sce", "ref_dolgalev_sce": "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_dolgalev_sce", "ref_lipka_sce": "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_lipka_sce"}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 14, "cluster": {}}
 cd /home/l012t/repositories/Interspecies_BM_phd/code/03_integration && \
/home/l012t/micromamba/envs/snakemake-tutorial/bin/python3.8 \
-m snakemake /omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_ctyp/sce_mspr-07 --snakefile /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/integration_snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.zb2uj81o' '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge/sce_mspr-06' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /home/l012t/repositories/Interspecies_BM_phd/code/config/config.yaml -p --allowed-rules cell_type_annotation --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /home/l012t/micromamba/envs/snakemake-tutorial/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.zb2uj81o/14.jobfinished || (touch /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.zb2uj81o/14.jobfailed; exit 1)

