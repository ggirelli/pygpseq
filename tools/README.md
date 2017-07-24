pygpseq-tools
=============

A collection of scripts that can be used in combination with pygpseq to perform advanced operations.

Each script needs to be able to load the `pygpseq` library. Either install it systemwide or place a link to the `pygpseq` folder in the same directory as the script (e.g., `ln -s REPOFOLDER/pygpseq ./pygpseq`).

* **analyze_dots.py** requires the table output of DOTTER, with the FISH dots coordinates in the `x`, `y` and `z` columns, and the NM address in the `File` column. Then, prepare a folder containing the DNA staining channels (deconvolved) with file name in DOTTER notation (e.g., `dapi_001_cmle.tif`). The script will add 5 columns: `cellID`, `lamin_dist`, `lamin_dist_norm`, `centr_dist` and `centr_dist_norm`. Also, an option is available to specify the voxel aspect ratio. Use `./analyze_dots.py` for more details.