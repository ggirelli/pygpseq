How to
======

Set up your data for pyGPSeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**pyGPSeq** requires your data (images) to be already in `.tif` format. The following steps include also convertion of your images to `.tif` format using the `toTif_fiji.js` script.

1. First of all, create an empty directory which will contain your data. We will refer to this diretory as **basedir**.
2. Your `.tif` images should be moved to subfolders of **basedir**. Note that each subfolder in **basedir** is expected to be a different condition. Thus, group in the same subfolder images from the same condition, and divide in different subfolders images from different conditions.
    * If your images are not in `.tif` format but in a microscope-specific format, you can use the `toTif_fiji.js` script. The script will convert your images to `.tif` format (one image per channel), moving them in a subfolder with the same name as the initial image.
    * Given how the script works, it is convenient to have one microscope-specific formatted image per condition. If that is the case, one subfolder will be created for every condition, containing the `.tif` images in a proper format. These subfolders can then be moved to **basedir** to proceed with the analysis.

Run pyGPSeq
~~~~~~~~~~~

**pyGPSeq** can be run on a single dataset or on multiple datasets (batch) at once. In the latter case, a queue is created and only one dataset at a time is analyzed.