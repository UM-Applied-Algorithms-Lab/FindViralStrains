# FindViralStrains
The goal of this pipeline is to sift out viral strains from data sets using MFD.
These are early drafts and ideas for our workflow, and are not yet polished and/or finished.
If you do have questions about anything that is happening in this workflow feel free to reach
out to us via github or email. Below you'll find the basic setup instructions that have been
tested on Redhat and Debian based Linux systems. 

Cheers,
McKayl & Lucy


Please note that everything below assumes that you already have copied our code into a working
directory, and have some variation of Conda installed on your device. If you need help doing any of
this below i've linked some very informative guides. I would suggest reaching out to your system
administrators for help/training if you are setting this up on any sort of HPC enviroment.

https://github.com/git-guides/git-clone
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

To set up your own working enviroment, I would use our CondaTemplate.yml file, along with the
command below.

```conda env create -f CondaTemplate.yml```

Once you have build list and verified that all packages are successfully installed you can move on to
actually running 

```snakemake -s findviralstrains_noref.smk --configfile config_files/no_ref_test.yml --cores 2```

Note: `cuttlefish` needs to be able to create many temporary files; you probably need to run
```
ulimit -n 2048
```
to allow your system to create so many files. If not, you will get this error in the cuttlefish step:
```
Error: Cannot open temporary file output/Test1/kmc_01021.bin:
```
