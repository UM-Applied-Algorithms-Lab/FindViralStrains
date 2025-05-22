# FindViralStrains
The goal of this pipeline is to sift out viral strains from data sets using MFD.
These are early drafts and ideas for our workflow, and are not yet polished and/or finished.
If you do have questions about anything that is happening in this workflow feel free to reach
out to us via github or email. Below you'll find the basic setup instructions that have been
tested on Redhat and Debian based Linux systems.

Cheers,
McKayl, Lucy, & Tim


Please note that everything below assumes that you already have copied our code into a working
directory, and have some variation of Conda installed on your device. If you need help doing any of
this below i've linked some very informative guides. I would suggest reaching out to your system
administrators for help/training if you are setting this up on any sort of HPC enviroment. The
command to download our repository is below.

```
git clone --recurse-submodules https://github.com/UM-Applied-Algorithms-Lab/FindViralStrains
```
https://github.com/git-guides/git-clone
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

Set up your own working enviroment, I would use our environment.yml file, along with the
command below. As written, this will create an environment called `find_viral_strains`---you can change this by editing
the first line of `environment.yml`.

```
conda env create -f environment.yml
```

Now, compile the rust libraries that we use.

```
cargo build --release
```

Once you have build list and verified that all packages are successfully installed you can move on to
actually running. An example config file is `config_files/example_config.yaml`.
If you have cloned the directory, it should run as is.

```
snakemake -s findviralstrains.smk --configfile config_files/example_config.yaml --cores 2
```

# How to contribute to this repository
Todo, but generally: make a branch for your new feature, make changes, commit them, merge into main, then make a pull request.
