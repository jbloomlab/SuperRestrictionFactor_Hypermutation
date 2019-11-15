# Analyzing mutation signitures of an APOBEC super restriction factor
### Adam S. Dingens, Jesse Bloom
### In collabortation Mollie McDonnell and Michael Emerman


FILL IN 

Computational analysis performed by Adam Dingens in the [Bloom lab](http://research.fhcrc.org/bloom/en.html) in Fall 2019. 


## Running the analysis
The main analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).

To run [analysis_notebook.ipynb](analysis_notebook.ipynb) and generate the [Markdown results](results/analysis_notebook.md), run the bash script [run_notebook.bash](run_notebook.bash) with:

./run_notebook.bash

On the Hutch cluster, you probably want to submit this job using [slurm](https://slurm.schedmd.com/), which you can simply do with:

sbatch -p largenode -c 16 --mem=100000 run_notebook.bash


## Organization
The mutational antigenic profiling analysis is performed by the iPython notebook [`analysis_notebook.ipynb`](analysis_notebook.ipynb). 

Subdirectories:

   * `./results/` iss generated by [`analysis_notebook.ipynb`](analysis_notebook.ipynb). 
   
   * `./results/FASTQ_files/` contains the input Illumina deep sequencing data. This file is generated by [`analysis_notebook.ipynb`](analysis_notebook.ipynb), which downloads the sequencing files from the [Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra).

   * `./data/` contains all input data needed to run [`analysis_notebook.ipynb`](analysis_notebook.ipynb)). Files are described in the iPython notebok when used. 


## Results and Conclusions
All of the results are breifly detailed in the [`analysis_notebook.ipynb`](analysis_notebook.ipynb) notebook; click on this notebook to see the results.



