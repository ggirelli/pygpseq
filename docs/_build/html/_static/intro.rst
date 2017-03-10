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

This step of the pipeline iterates over the conditions detected in the beginning. For every condition, the pixels of every nucleus are retrieved, specifically, the following pixel information are obtained: series and nucleus id number, DNA channel intensity (minus background), signal channel intensity  (minus background), distance from nuclear lamina and relative distance from nuclear lamina.

4. General boxplots
~~~~~~~~~~~~~~~~~~~

5. Final plots
~~~~~~~~~~~~~~

6. Final report
~~~~~~~~~~~~~~~

References
----------

.. [1] Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." Automatica 11.285-296 (1975): 23-27.
.. [2] http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_adaptive