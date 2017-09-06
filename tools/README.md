pygpseq-tools
=============

A collection of scripts that can be used in combination with pygpseq to perform advanced operations.

Each script needs to be able to load the `pygpseq` library. Either install it systemwide or place a link to the `pygpseq` folder in the same directory as the script (e.g., `ln -s REPOFOLDER/pygpseq ./pygpseq`).

* **dotter2gpseq.py** requires the table output of DOTTER, with the FISH dots coordinates in the `x`, `y` and `z` columns, and the NM address in the `File` column. Then, prepare a folder containing the DNA staining channels (deconvolved) with file name in DOTTER notation (e.g., `dapi_001_cmle.tif`). The script will add 5 columns: `cellID`, `lamin_dist`, `lamin_dist_norm`, `centr_dist` and `centr_dist_norm`. Also, an option is available to specify the voxel aspect ratio. Moreover, if you use the `-s` option, the script identifies each nucleus in the images and selects those in G1 based on the distribution of their flattened size and DNA stain integral. The G1 selection works only if the images are enough populated (at least ~300 nuclei) and is relatively time consuming. Use `./analyze_dots.py` for more details.