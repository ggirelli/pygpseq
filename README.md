pyGPSeq v2.0.0
=======================

A Python package that provides tools to analyze images of GPSeq samples.

Sample scripts are available showing the code for single/batch runs.

Read the [documentation](https://github.com/ggirelli/gpseq-img-py/wiki) for more details.

Installation
-------------

```
git clone http://github.com/ggirelli/gpseq-seq-py
cd gpseq-seq-py
sudo -H pip3 install -e .
```

Usage
----------

Use the `gpseq_anim` script.

```
usage: gpseq_anim [-h]
                  [--skip [{inst,seg,an,box,plot,report} [{inst,seg,an,box,plot,report} ...]]]
                  [-l log] [-a Z Y X] [-U unit] [-d dna_name [dna_name ...]]
                  [-s sig_name [sig_name ...]] [-z min_z]
                  [--seg-type {sum_proj,max_proj,3d}]
                  [--an-type {sum_proj,max_proj,3d,mid}]
                  [--mid-type {central,largest,maxIsum}]
                  [--nuclear-sel [{size,surf,shape,sumI,meanI,flat_size} [{size,surf,shape,sumI,meanI,flat_size} ...]]]
                  [--sigma-smooth sigmaValue] [--sigma-density sigmaValue]
                  [--description [DESCRIPTION [DESCRIPTION ...]]]
                  [-t nthreads] [--note NOTE] [--regexp REGEXP] [-r] [-n] [-u]
                  [--version]
                  inDir outDir

Run GPSeq image analysis.

positional arguments:
  inDir                 Path to input directory, containing singlcde-condition
                        directories with TIF files.
  outDir                Path to output directory, must be different from the
                        input directory.

optional arguments:
  -h, --help            show this help message and exit
  --skip [{inst,seg,an,box,plot,report} [{inst,seg,an,box,plot,report} ...]]
                        Space-separated phases to be skipped. Use -- after the
                        last one.
  -l log, --logpath log
                        Path to log file. By default: outDir/log
  -a Z Y X, --aspect Z Y X
                        Physical size of Z, Y and X voxel sides. Default:
                        300.0 216.6 216.6
  -U unit, --umes unit  Unit of measure for the aspect. Default: nm
  -d dna_name [dna_name ...], --dna-channels dna_name [dna_name ...]
                        Space-separated names of DNA staining channels. Use --
                        after the last one.
  -s sig_name [sig_name ...], --sig-channels sig_name [sig_name ...]
                        Space-separated names of GPSeq signal channels. Use --
                        after the last one.
  -z min_z, --min-z min_z
                        If lower than 1, minimum fraction of stack, if higher
                        than 1, minimum number of slices to be occupied by a
                        nucleus. Default: .25
  --seg-type {sum_proj,max_proj,3d}
                        Segmentation type. Default: 3d
  --an-type {sum_proj,max_proj,3d,mid}
                        Analysis type. Default: mid
  --mid-type {central,largest,maxIsum}
                        Method for mid-section selection. Default: largest
  --nuclear-sel [{size,surf,shape,sumI,meanI,flat_size} [{size,surf,shape,sumI,meanI,flat_size} ...]]
                        Space-separated features for nuclear selection. Use --
                        after the last one. Default: flat_size sumI
  --sigma-smooth sigmaValue
                        Sigma value for sparse gaussian smoothing.
  --sigma-density sigmaValue
                        Sigma value for density calculation.
  --description [DESCRIPTION [DESCRIPTION ...]]
                        Space separated condition:description couples.
                        'condition' are the name of condition folders.
                        'description' are descriptive labels used in plots
                        instead of folder names. Use -- after the last one.
  -t nthreads, --threads nthreads
                        Number of threads to be used for parallelization.
                        Increasing the number of threads might increase the
                        required amount of RAM.
  --note NOTE           Dataset/Analysis description. Use double quotes.
  --regexp REGEXP       Advanced. Regular expression to identify tif images.
  -r, --rescale-deconvolved
                        Perform rescaling of deconvolved images. Requires
                        Huygens Professional v4.5 log file for an image to be
                        rescaled.
  -n, --normalize-distance
                        Perform distance normalization. Necessary to compare
                        nuclei with different radius.
  -u, --DEBUG-MODE      Debugging mode.
  --version             show program's version number and exit
```

License
---

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```