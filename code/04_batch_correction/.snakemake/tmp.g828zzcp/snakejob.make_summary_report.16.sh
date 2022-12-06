#!/bin/sh
# properties = {"type": "single", "rule": "make_summary_report", "local": false, "input": ["/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge", "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_rnrm", "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/08_mnncorrect"], "output": ["/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/reports/integration/integration_summary.html"], "wildcards": {}, "params": {"nr_hvgs": 2000, "species_all": ["mcar", "mcas", "mmus", "mspr", "all"], "batches": ["Batch_exp_day"]}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 16, "cluster": {}}
 cd /home/l012t/repositories/Interspecies_BM_phd/code/03_integration && \
/home/l012t/micromamba/envs/snakemake-tutorial/bin/python3.8 \
-m snakemake /omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/reports/integration/integration_summary.html --snakefile /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/integration_snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.g828zzcp' '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge' '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_rnrm' '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/08_mnncorrect' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /home/l012t/repositories/Interspecies_BM_phd/code/config/config.yaml --config run_integration_summary=True -p --allowed-rules make_summary_report --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /home/l012t/micromamba/envs/snakemake-tutorial/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.g828zzcp/16.jobfinished || (touch /home/l012t/repositories/Interspecies_BM_phd/code/03_integration/.snakemake/tmp.g828zzcp/16.jobfailed; exit 1)

