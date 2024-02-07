# RNASeq-DE-nf

<p align="center">
:wrench: This pipeline is currently under development and is not currently functional :wrench:
</p> 

  - [Description](#description)
  - [Diagram](#diagram)
  - [How to cite this workflow](#how-to-cite-this-workflow)
  - [User guide](#user-guide)
      - [Quick start guide](#quick-start-guide)
      - [Install instructions](#install)
      - [Dependencies & third party tools](#dependencies--third-party-tools)
      - [Required (minimum)
        inputs/parameters](#required-minimum-inputsparameters)
      - [Recommendations for use on specific compute systems](#recommendations-for-use-on-specific-compute-systems)
      - [Compute resource usage on tested
        infrastructures](#compute-resource-usage-across-tested-infrastructures)
      - [Benchmarking (compute resource usage on tested infrastructures)](#benchmarking--compute-resource-usage-on-tested-infrastructures-)
  - [Additional notes](#additional-notes)
  - [Help/FAQ/Troubleshooting](#helpfaqtroubleshooting)
  - [3rd party Tutorials](#3rd-party-tutorials)
  - [Licence(s)](#licences)
  - [Acknowledgements/citations/credits](#acknowledgementscitationscredits)

---

## Description

RNASeq-DE-nf is a pipeline that pre-processes Illumina RNA sequencing data for differential expression (raw FASTQ to counts). The pipeline is written in Nextflow and uses Singularity to run containerised tools.  

The steps in this pipeline include:
1. QC of raw FASTQs (Fastqc -> multiQC)
2. Trim raw FASTQs (bbduk)
3. Mapping with STAR aligner 
  -  Alignment (Fastq-to-bam)
  -  Merge lane level to sample level BAMs
  -  Mapping metrics
  -  Raw counts: HTSeq (bam-to-counts)  
4. Salmon pseudo-alinement (Fastq-to-counts)

## Diagram

Logical visual description of processing steps for workflow (TBD)


## How to cite this workflow

Add citation instructions here.


## User guide

### 1. Prepare inputs
The scripts in this repository use relative paths and require careful setup to ensure that scripts can locate input files seamlessly. To start:

To run this pipeline, you will need to prepare your input files, reference data, and clone this repository. Before proceeding, ensure Nextflow is installed on the system you're working on. To install Nextflow, see these [instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation).  

To run this pipeline you will need the following inputs: 
- **Illumina raw FASTQ files**
- **Reference files**
    - Reference genome primary assembly (.fasta) and corresponding annotation (.gtf) file 
    - Annotation in BED format is required for RSeQC's infer_experiment.py (This file, if not provided can also be created using the `.gtf` file)
    - References can be obtained: 
      * following recommendations in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). 
      * from [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html)
      * Indices can also be Generated when excecuting the pipeline

- Input sample sheet (`.csv`) file

You will need to create a sample sheet with information about the samples you are processing, before running the pipeline. This file must be **comma-separated** and contain a header and one row per sample. Columns should correspond to sampleID, BAM-N file-path, BAM-T file-path: 

```csv
sampleID,Lane,R1,R2,SEQUENCING_CENTRE,PLATFORM,RUN_TYPE_SINGLE_PAIRED
SAMPLEID1,sample1_1.fastq.gz,sample1_2.fastq.gz,KCCG,ILLUMINA,PAIRED
SAMPLEID2,sample2_1.fastq.gz,sample2_2.fastq.gz,KCCG,ILLUMINA,PAIRED
SAMPLEID2,sample2.fastq.gz,,KCCG,ILLUMINA,SINGLE
```
When you run the pipeline, you will use the mandatory `--input` parameter to specify the location and name of the input sample sheet file: 
```
--input /path/to/samples.csv
```

### 2. Prepare the reference materials 



### 3. Clone this repository 

Download the code contained in this repository with: 

```
git clone https://github.com/Sydney-Informatics-Hub/RNAseq-DE-nf
```

This will create a directory with the following structure: 
```
Somatic-shortV-nf/
├── LICENSE
├── README.md
├── config/
├── images/
├── scripts/
├── main.nf
├── modules/
└── nextflow.config
```
The important features are: 

* **main.nf** contains the main nextflow script that calls all the processes in the workflow.
* **nextflow.config** contains default parameters to use in the pipeline.
* **modules** contains individual process files for each step in the workflow. 
* **config** contains infrastructure-specific config files (this is currently under development)


### 4. Run the pipeline 

The minimal run command for executing this pipeline is: 

```
nextflow run main.nf --input samples.csv \
                     --refFasta /path/to/ref.fasta --refGtf /path/to/ref.gtf
                     
```

By default, this will generate `work` directory, `results` output directory and a `runInfo` run metrics directory inside the results directory. 

To specify additional optional tool-specific parameters, see what flags are supported by running:

```
nextflow run main.nf --help 
```

The  nextflow command with optional parameters, where you can define the name of the output folder and provide an integer value for the number of genomic-interval files is:
```
To be updated
```
**Mandatory parameters** 

**Optional parameters**  


If for any reason your workflow fails, you are able to resume the workflow from the last successful process with `-resume`. 

This pipeline has been optimised for NCI Gadi HPC, instructions for executing on Gadi are provided in [Infrastructure usage and recommendations](#infrastructure-usage-and-recommendations)


### 5. Results 

## Infrastructure usage and recommendations
This pipeline has been successfully implemented on [NCI Gadi HPC](https://nci.org.au/our-systems/hpc-systems) using a infrastructure-specific config. 
As per the `config/gadi.config`, the main script is excecuted using the queue `copyq` so that the singularity container images required by the pipeline are downloaded. The NCI Gadi config currently runs all other tasks (except downloading the singularity images) on the normal queue. This config can be used to interact with the job scheduler and assign a project code to all task job submissions. 
The following flags are required to be specified in the command:

* `--whoami` your NCI or Pawsey user name
* `--gadi-account` the Gadi project account you would like to bill service units to

The config uses the `--gadi-account` flag to assign a project code to all task job submissions for billing purposes. The version of Nextflow installed on Gadi has been modified to make it easier to specify resource options for jobs submitted to the cluster. See NCI's [Gadi user guide](https://opus.nci.org.au/display/DAE/Nextflow) for more details.

The minimal run command for executing this pipeline on NCI Gadi HPC is:

```
nextflow run main.nf -resume \
        --input samples.csv \
        -profile gadi \
        --whoami $(whoami) --gadi_account $PROJECT \
        --NCPUS ${NCPUS} \
       --refFasta /path/to/ref.fasta --refGtf /path/to/ref.gtf
```

Before running the pipeline you will need to load Nextflow and Singularity, both of which are globally installed modules on Gadi. You can do this by running the commands below:

```
module purge
module load nextflow singularity
```




### Quick start guide

General guide for deployment across multiple infrastructures (distinct from specific infrastructure quick start guide)


### Install

General installation guide.

> If there are different installation requirements for specific compute infrastructures you could indicate these here, or in an individual infrastructure documentation template: https://github.com/AustralianBioCommons/doc_guidelines/blob/master/infrastructure_optimisation.md


### Dependencies & third party tools


### Required (minimum) inputs/parameters

he minimum inputs required for the workflow to run.


### Recommendations for use on specific compute systems

link to installation instructions for each infrastructure 
recommendations
    
Documentation for a specific infrastructure should go into a infrastructure documentation template
https://github.com/AustralianBioCommons/doc_guidelines/blob/master/infrastructure_optimisation.md


### Benchmarking
Coming soon!


## Additional notes


## Help / FAQ / Troubleshooting


## 3rd party Tutorials 


## [License(s)](../LICENSE.md)
GNU General Public License v3.0


## Acknowledgements/citations/credits

### Authors 
- Tracy Chew
- Rosemarie Sadsad
- Ching-Yu Lu (Sydney Informatics Hub, University of Sydney)   
- Nandan Deshpande (Sydney Informatics Hub, University of Sydney)
- Georgie Samaha (Sydney Informatics Hub, University of Sydney)

### Acknowledgements 

- This pipeline was built using the [Nextflow DSL2 template](https://github.com/Sydney-Informatics-Hub/Nextflow_DSL2_template).  
- Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).  

### Cite us to support us! 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:  

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia. 
