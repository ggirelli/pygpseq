pygpseq-tools
=============

A collection of scripts that can be used in combination with pygpseq to perform advanced operations.

Each script needs to be able to load the `pygpseq` library. Either install it systemwide or place a link to the `pygpseq` folder in the same directory as the script (e.g., `ln -s REPOFOLDER/pygpseq ./pygpseq`).

### dotter2gpseq.py

**dotter2gpseq.py** requires the table output of DOTTER, with the FISH dots coordinates in the `x`, `y` and `z` columns, and the NM address in the `File` column. Then, prepare a folder containing the DNA staining channels (deconvolved) with file name in DOTTER notation (e.g., `dapi_001_cmle.tif`).

The script will add 7 columns: `angle`, `Allele`, `cellID`, `lamin_dist`, `lamin_dist_norm`, `centr_dist` and `centr_dist_norm`.

The `Allele` column contains allele labeling:

- No value: dot outside of cells
- -1: more than 2 dots per cell/channel
- 0: less than 2 dots per cell/channel
- 1: more central dot
- 2: more peripheral dot

The `angle` column contains the angle between each allele pair and the nucleus center of mass.

Lamina distance and center distance are calculated using an anisotropic euclidean transform based on the specified aspect (`-a`). Lamina is defined as the first black (background) pixels outside of the nucleus. Center is defined as the nuclear pixels most distant from the lamina. Previous versions of the script extrapolated the center distance from the normalized lamin distance instead.

Distances from center and lamina are normalized on their sum, for each dot. Previous versions of the script normalized over the maximum radius of the nucleus, *assuming a spherical shape*.

* Also, an option (`-a Z Y X`) is available to specify the voxel aspect ratio.
* Use `-t` to specify the number of threads to be used for paralelization.
* Use `--dilate` to specify the number of dilation operations to perform. The dilations are performed nucleus-wise, while the nuclear masked is saved for a general dilation (for simplicity).
* Use `./dotter2gpseq.py -h` for more details.