<h1 align = 'center'>🧬Rachel's <code>Bioinformatics Pipelines</code></h1>
<h2 align = 'center'>⚠️A word of warning...</h2>
This repository contains 'omics data-type based pipelines from my general role as a bioinformatician. These are not cleaned, nor are they transportable. However, I figured they are worth allowing others to use if they want!

<h2 align = 'center'>▶️ Usage</h2>

- To run many of these scripts, users must have access to a high performance computing environment. They are often in pairs with the same name, with the shell script (ex. shell.sh) created to submit the working script (shell.R) to an HPC environment.
- The user must update the bash script to the appropriate parameters for their unique HPC. These were built on an HPC using command line with SLURM.

<h2 align = 'center'>🚀Current pipelines</h2>

| Pipeline                              | Description                                                                                                                                                                                                 |
|---------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 🧬 `Bulk RNAseq`                      | Basic pipeline for bulk RNA-sequencing experiments from preprocessing to downstream pathway analysis                                                                                                        |
| 📦 `PackageBuilding`                  | Helper scripts for building an R package                                                                                                                                                                    |
| 📊 `Visualizations`                   | Basic visualization scripts for common plots such as volcano plots                                                                                                                                          |
| 🎯 `miRNA`                            | Scripts for analysis of miRNA                                                                                                                                                                               |
| 🧬 `scDNAseq/MissionBio Tapestri`     | Analysis scripts for MissionBio Tapestri platform data using the `optima` package and hardcoding analysis                                                                                                   |
| 🔍 `scRNAseq`                         | Basic analysis for scRNAseq from 10X platforms                                                                                                                                                              |
| 🧭 `spatial RNAseq/VisiumHD`          | Scripts from analyses of Visium HD samples for preprocessing, quality control, integration using `scVI`, nucleus segmentation using `stardist`, differential expression analysis, cell type ID via RCTD    |


---

<h2 align = 'center'>🗓️ Last Updated</h2>

📅 **June 24, 2025**  
✍️ Maintained by [Rachel Griffard-Smith](https://github.com/rachelgriffard)  
📫 Contact: Leave a note in the `Issues` or `Pull requests` section
