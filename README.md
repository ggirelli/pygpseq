pyGPSeq
=======================

A Python package that provides tools to analyze images of GPSeq samples.

Sample script are available showing the code for single runs and batch runs.

Read the documentation in ``docs/_build/html`` for more details.

Build package
-------------

Source distribution: `python setup.py sdist`
Universal wheel: `python setup.py bdist_wheel --universal`

Build docs
----------

From withing the root folder: `sphinx-apidoc -f -o docs/ .`

How to run
----------

Use the `gpseq_anim` script.

```
usage: gpseq_anim  [-h] [--skip {inst,seg,an,box,plot,report}] [-l log]
                   [-a Z Y X] [-d [dna_name [dna_name ...]]]
                   [-s [sig_name [sig_name ...]]] [-z min_z]
                   [--seg-type {sum_proj,max_proj,3d}]
                   [--an-type {sum_proj,an_proj,3d,mid}]
                   [--mid-type {central,largest,maxIsum}]
                   [--nuclear-sel [{size,surf,shape,sumI,meanI,flat_size} [{size,surf,shape,sumI,meanI,flat_size} ...]]]
                   [--description [DESCRIPTION [DESCRIPTION ...]]]
                   [-t ncores] [--note NOTE] [--regexp REGEXP] [-r] [-n]
                   inDir outDir

Run GPSeq image analysis.

positional arguments:
  inDir                 Path to input directory, containing single-condition
                        directories with TIF files.
  outDir                Path to output directory, must be different from the
                        input directory.

optional arguments:
  -h, --help            show this help message and exit
  --skip {inst,seg,an,box,plot,report}
                        Space-separated phases to be skipped. Use -- after the
                        last one.
  -l log, --logpath log
                        Path to log file. By default: outDir/log
  -a Z Y X, --aspect Z Y X
                        Physical size of Z, Y and X voxel sides. Default:
                        300.0 216.6 216.6
  -d [dna_name [dna_name ...]], --dna-channels [dna_name [dna_name ...]]
                        Space-separated names of DNA staining channels. Use --
                        after the last one.
  -s [sig_name [sig_name ...]], --sig-channels [sig_name [sig_name ...]]
                        Space-separated names of GPSeq signal channels. Use --
                        after the last one.
  -z min_z, --min-z min_z
                        If lower than 1, minimum fraction of stack, if higher
                        than 1, minimum number of slicesto be occupied by a
                        nucleus
  --seg-type {sum_proj,max_proj,3d}
                        Segmentation type. Default: 3d
  --an-type {sum_proj,an_proj,3d,mid}
                        Analysis type. Default: mid
  --mid-type {central,largest,maxIsum}
                        Method for mid-section selection.
  --nuclear-sel [{size,surf,shape,sumI,meanI,flat_size} [{size,surf,shape,sumI,meanI,flat_size} ...]]
                        Space-separated features for nuclear selection. Use --
                        after the last one. Default: flat_size sumI
  --description [DESCRIPTION [DESCRIPTION ...]]
                        Space separated condition:description couples.
                        'condition' are the name of condition folders.
                        'description' are descriptive labels used in plots
                        instead of folder names. Use -- after the last one.
  -t ncores, --threads ncores
                        Number of threads to be used for parallelization.
                        Increasing the number of threads might increase the
                        required amount of RAM.
  --note NOTE           Dataset/Analysis description.
  --regexp REGEXP       Advanced. Regular expression to identify tif images.
  -r, --rescale-deconvolved
                        Perform rescaling of deconvolved images. Requires
                        Huygens Professional v4.5 log file for an image to be
                        rescaled.
  -n, --normalize-distance
                        Perform distance normalization. Necessary to compare
                        nuclei with different radius.
```
