## Dependencies:
* Ubuntu 20.04 or WSL2
* Python 2.7 
* Dependency packages such as bwa, bedtools, numpy, matplotlib (in environment.yml)
* Reference genome fasta file(We have uploaded the reference genome fasta to Baidu web disk. [human](https://pan.baidu.com/s/1wKNr_JJuuuvWkngL_jtiqQ):   cpfa &nbsp;&nbsp;&nbsp;[mouce](https://pan.baidu.com/s/1UxwSA0Rx_PqdvaX-kNPqFg):z1ce)
## Installation:
    # We recommend doing the following in a conda. A new conda environment will be created with Environment.yml.
    # Assume that you already own conda and have a conda base environment open...
    wget https://github.com/qwe1234567891/OliTag-seq/archive/refs/heads/main.zip
    unzip main.zip  # Download the olitag code package and unzip it...
    cd OliTag-seq-main
    
    conda env create -f environment.yml  # Create a new conda environment for olitag...
    source activate olitag
## Example and description of the Manifest File(\*.yaml)
We added the input parameters required by olitag, the reference genome path, and some tool information to manifest.yaml so that we could manage and adjust the input parameters manifest.yaml can be opened and edited using most text editing tools. Here is an example:

    reference_genome: ../reference_genome/human.fa
    output_folder: ./

    bwa: bwa
    bedtools: bedtools

    data1: ["test_1.fq.gz"]
    data2: ["test_2.fq.gz"]


    samples:
      EXM1:
        target: GTGGGGAAGAAGGTGTCTTCNGG
        barcode1: AACCTCTT
        barcode2: CCAATCTG
        description: EXM1

      EXM2:
        target: GTGGGGAAGAAGGTGTCTTCNGG
        barcode1: AATACCGC
        barcode2: CCAATCTG
        description: EXM2
Meaning of each field:
* _reference_genome:_ Reference to the genomic fasta path can be human or mouse.
* _output_folder:_ Olitage-seq Specifies the path for storing generated files.
* _bwa、bedtools:_ The sequence alignment tool used in the code.
* _data1、data2:_ The input data(\*_1.fq.gz、\*_2.fq.gz)
* _samples:_ A nested field containing the details of each sample. The required parameters are targetsites, forward barcode reverse barcode and a descroption of the sample.

## Running:
After downloading our code package and installing it successfully, you can run our test data using the following commands:

    python Lentegrate-Seq-mouse/Lentegrate-Seq.py all -m manifest.yaml
