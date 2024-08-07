# FindViralStrains
The goal of this pipeline is to sift out viral strains from data sets using MFD.
These are early drafts and ideas for our workflow, and are not yet polished and/or finished.
If you do have questions about anything that is happening in this workflow feel free to reach
out to us via github or email. Below you'll find the basic setup instructions for a Redhat
based Linux system.

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

>conda env create -f CondaTemplate.yml

Once you have build list and verified that all packages are successfully installed you can move on to
actually running 
