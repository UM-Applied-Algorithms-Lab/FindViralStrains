# FindViralStrains
The goal of this pipeline is to sift out viral strains from data sets using MFD.
Some of these programs may not yet polished and/or completely finished.
If you do have questions about anything that is happening in this workflow feel free to reach
out to us via github or email. Below you'll find the basic setup instructions that have been
tested on Redhat, Debian, and Mac OS.

Cheers,
University of Montana Applied Algorithms Lab

First, you should clone our repository to your local machine using the command below.

```
git clone --recurse-submodules https://github.com/UM-Applied-Algorithms-Lab/FindViralStrains
```

Next you will need to download dependencies and set up your local enviroment. We will be
using Conda/Anaconda for this. If you do not have this locally you can use the guides
below. 

https://github.com/git-guides/git-clone
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

This command downloads program packages from our pre-configured file that will be needed 
to run the pipeline. 

```
conda env create -f environment.yml
```

Now, compile the rust libraries that we use.

```
cargo build --release
```

Additionally, you will need to locally donwload the emboss package using brew, apt, or
other similar package managers. Below is documentation for Mac users who do not already
have brew.

https://docs.brew.sh/Installation

```
brew install brewsci/bio/emboss
```

With all of the needed software now download, all that is left it so configure your config
file, and run the pipeline. In the config_files directory, you can edit the example file
to point to your data. Using the command below you can now run the entire pipeline.

```
snakemake -s findviralstrains.smk --configfile config_files/example_config.yaml --cores 2
```

After running through all the Snakemake rules, your output will be organized into the following directories:

```
decomp_results (paths chosen by the ILP solver and their visualizations)
graphs (data files for the graphs themselves)
output_genomes (final genome from each found path)
read_data (intermittent data files from your fastq files)
```

# How to contribute to this repository
Todo, but generally: make a branch for your new feature, make changes, commit them, merge into main, then make a pull request.
