Introduction
=====

**pygpseq** is a Python (**py**) package that provides the tools to analyze microscopy images of **GPSeq** samples. **GPSeq** (Genomic loci Positioning by Sequencing) is a protocol developed in the Bienko lab and, as per title, allows to estimate the radial positioning of genomic loci through sequencing. The imaging step, prior to sequencing, represents an informative step on the quality of the samples and provides additional information for the sequencing data analysis.

Pipeline steps
--------------

1. **Initiation**: identify image series and conditions.
2. **Segmentation**: identify nuclei in the images.
3. **Single-nucleus analysis**: retrieve nuclear pixels and calculate profiles.
4. **General boxplots**: study per-pixel and per-nucleus feature distributions.
5. **Final plots**: compare different conditions.
6. **Final report**: produce general pdf report.

1. Initiation
~~~~~~~~~~~~~

The pipeline starts by looking into the specified ```basedir``` (or base directory) for the GPSeq data. The data is expected to be divided in multiple sub-directories, each of which contains images of samples treated with different conditions. Thus, every sub-directory in `basdir` is considered to be a **condition directory**.
::

	- basedir
	  + condition #1
	    * series1.channel1.tif
	    * series1.channel2.tif
	    * series2.channel1.tif
	    * series2.channel2.tif
	  + condition #2
	    * ...
	  + ...


The extension of the images can be specified in the ```ext``` attribute of the class ```Main```, the default is ```.tif```.

The format of the image name can be specified as a regular expression in the ```reg``` attribute of the calls ```Main```. The default is reported in the ```Main``` class documentation.

In other words, the first step consists in the identification of the condition directories containing properly named tif images, to then proceed with their analysis.

2. Segmentation
~~~~~~~~~~~~~~~

The pipeline continues by going through the series in every condition and identifying the nuclei. Moreover, the median intensity of the background is calculated to be used later during the background removal in the single-nucleus analysis step (3). The segmentation is always performed on the channel containing the DNA staining, and the generated mask is used to estimate the background for the channel containing the linker signal.

The segmentation is performed as follows:

* Identify a (stack) global intensity threshold value using the Otsu's method. [1]_
* Binarize the image using the global intensity threshold value.
* Binarize the image using an adaptive (local) intensity threshold algorithm that compares the median value of the neighbourhood of every pixel with its intensity. [2]_ The size of the side of the square neighbourhood can be specified in the ```adaptive_neighbourhood``` attribute of the ```Main``` class.
* Combine the global and local masks (logical ```AND``` operation).

Then, the masked objects that touch the image XY borders are removed, and holes are filled both in 3D and on every XY plane (i.e., slice).

Finally, the mask is passed through some filters:

* A first filter on the size in the XY plane is applied to remove small objects from the mask.
* If performing a 3D analysis, the masked objects are filtered for their size on the Z dimension.

Before proceeding with the single nucleus analysis, the objects (i.e., nuclei) in the masked are stored in the ```Nuclei``` class. Specifically, the following information are retained:

* Series and nucleus id number (1-indexed).
* The bounding square (2D) or box (3D).
* The pixel/voxel sides ratio.
* The background in the original image channels.
* The size both in 2D and in 3D (if available).
* The surface (only in 3D and if turned on in the ```Main``` class).
* The sum and mean of the pixel intensities.
* The shape of the nucleus (e.g., circularity, sphericity).
* The applied global threshold.

3. Single-nucleus analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

This step of the pipeline iterates over the conditions detected in the beginning. For every condition, the nuclei are selected based on the distribution of user-selected features (based on the FWHM of the highest peak) and the pixels of every selected nucleus are retrieved.

The features that can be selected for the nuclear selection are:

* Size: area if 2D and volume if 3D.
* Surface: available only for 3D segmentation.
* Shape: circularity if 2D, volume if 3D.
* Intensity sum: integral of DNA stain intensity.
* Intensity mean: average DNA stain intensity.
* Flattened size: area in the Z projection.

Specifically, the following pixel information are obtained: series and nucleus id number, DNA channel intensity (minus background), signal channel intensity  (minus background), distance from nuclear lamina and relative distance from nuclear lamina.

4. General boxplots
~~~~~~~~~~~~~~~~~~~

General boxplot are generated to recap different distributions for each condition. The following boxplots are generated:

* Cell size: area if 2D, volume if 3D.
* Cell shape: circularity if 2D, sphericity if 3D.
* Cell average DNA stain intensity.
* Cell sum (integral) of DNA stain intensity.
* Single-pixel DNA stain intensity.
* Single-pixel GPSeq signal intensity.

5. Final plots
~~~~~~~~~~~~~~

Then, the following additional plots are generated:

* A *pixel study* containing the distribution of the single-pixel values is generated for the two channels and their ratio. The study contains 8 plots, from top-to-bottom, left-to-right:
    - Pixel intensity distribution (boxplot) per binned relative distance from nuclear lamina.
    - Mean pixel intensity (both raw and smoothened) against relative distance from nuclear lamina.
    - Median pixel intensity (both raw and smoothened) against relative distance from nuclear lamina.
    - Mode pixel intensity (both raw and smoothened) against relative distance from nuclear lamina.
    - Number of pixels per relative distance from nuclear lamina.
    - Chalk-on-Blackboard plot. Which is a heatmap showing the binned distribution of channel intensity against a binned relative distance from nuclear lamina. Also, the smoothened mean, median and mode curves are reported.
    - Standard deviation of pixel intensity (both raw and smoothened) against relative distance from nuclear lamina.
* A *profiles* plot containing six plots, from left-to-right, top-to-bottom:
    - The DNA staining mean, median and mode smoothened curves.
    - The DNA staining mean, median and mode smoothened curves first derivative.
    - The GPSeq signal mean, median and mode smoothened curves.
    - The GPSeq signal mean, median and mode smoothened curves first derivative.
    - The channel ratio mean, median and mode smoothened curves.
    - The channel ratio mean, median and mode smoothened curves first derivative.

6. Final report
~~~~~~~~~~~~~~~

A single pdf file report is generated containing the input parameters and most of the generated plots.

References
----------

.. [1] Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." Automatica 11.285-296 (1975): 23-27.
.. [2] http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_adaptive