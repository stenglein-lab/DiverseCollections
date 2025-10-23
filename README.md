# DiverseCollections
A repository to house the scripts, appropriate data, and other associated content for the paper titled: Diverse collection of wild Drosophila melanogaster reveals the evolutionary history of a persistent viral infection


Here is an overview of how scripts and data files in this repository map to figures and tables in the paper:

- **Figure 1**: Maps were created using [create_sample_map.R](analyses/scripts/create_sample_map.R) using data files referenced therein.  In general, scripts refer to file paths relative to the root of this repository.  This can be created by setting the working directory to Project Directory in RStudio.  qPCR data panels were created using [qPCR_data.R](analyses/scripts/qPCR_data.R) using input datafiles referenced therein.  The map images and qPCR data panels were combined using [Affinity Designer](figures/affinity_designer/Figure_1.afdesign) and output as [a PDF](figures/Figure_1.pdf).

- **Supplemental Table 1**: is just a list of primer sequences.

- **Supplemental Tables 2 and 3**: were output from the [qPCR_data.R](analyses/scripts/qPCR_data.R) script.

- **Read preprocessing** was performed using our lab's [read_preprocessing pipeline](https://github.com/stenglein-lab/read_preprocessing) [[v1.0](https://github.com/stenglein-lab/read_preprocessing/releases/tag/v1.0)].  This pipeline trims adapter and low quality sequences from the ends of Illumina reads. We ran this using the command: `nextflow run stenglein-lab/read_preprocessing -profile singularity -resume  --fastq_dir collected_fastq --outdir preprocessed_reads`.

- **Supplemental Figure 1**: coverage plots were output using our lab's remapping_workflow, a nextflow workflow to map reads to a set of reference sequences.  The version of this workflow that was used for this paper is archived in the [virus_evol_paper branch of this repository](https://github.com/stenglein-lab/remapping_workflow/blob/virus_evol_paper/).  The entry point for this workflow was [this script](https://github.com/stenglein-lab/remapping_workflow/blob/virus_evol_paper/run_remapping_workflow).  The coverage plot figure was created using this [R script](https://github.com/stenglein-lab/remapping_workflow/blob/virus_evol_paper/bin/plot_refseq_coverage.R)

- **Supplemental Table 4**: is a table of metadata for sequenced samples

- **Supplemental Table 5**: is a table of species assignments based on drosophilid Cytochrome C oxidase subunit 1 sequences.  These assignments were generated using our lab's [species_id pipeline](https://github.com/stenglein-lab/species_id) [[v1.0](https://github.com/stenglein-lab/species_id/releases/tag/v1.0)], using the command line: `nextflow run stenglein-lab/species_id -profile singularity -resume  --fastq_dir ../preprocessed_reads/trimmed_fastq --outdir ../species_id -r v1.0 --fastq_pattern "*_[12].fastq.gz"

- 



