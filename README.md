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

To update, run the following from within the repository folder.

```
git pull
sudo -H pip3 install -e .
```

Usage
----------

#### Analyze a GPSeq image dataset

The `gpseq_anim` (**GPSeq** **an**alysis of **im**ages) analyzes a multi-condition GPSeq image dataset. Run `gpseq_anim -h` for more details.

#### Calculate lamin distance of FISH signals

The `gpseq_fromfish` script characterizes FISH signals identified with `DOTTER` (or similar tools) by calculating: absolute/normalized distance from lamina and central region, nuclear compartment, allele status,... Run `gpseq_fromfish -h` for more details.

#### Perform automatic 3D nuclei segmentation

Run `tiff_auto3dseg -h` for more details on how to produce binary/labeled (compressed) masks of your nuclei staining channels

#### Identify out of focus (OOF) fields of view

Run `tiff_findoof -h` for more details on how to quickly identify out of focus fields of view. Also, the `tiff_plotoof` script (in R, requires `argparser` and `ggplot2`) can be used to produce an informative plot with the signal location over the Z stack.

#### Uncompress a tiff

To uncompress a set of tiff, use the `tiff_uncompress` command (`-h` for more details).

License
---

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```