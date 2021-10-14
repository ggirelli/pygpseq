# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased

## [3.4.1] - 2021-11-14
This release fixes dependency conflicts that were caused by an improper dependency setup.
As such, no functionality is changed. The same can be said for the way a user interacts with the package.
The main change is related to the new dependency manager: Poetry.

We would still recommend using our new [`radiantkit`](https://github.com/ggirelli/radiantkit) package instead.
Use `pygpseq==3.3.6` only if strictly necessary and at your own risk, as this package is not currently maintained.

### Changed
- Now supporting specifically only Python 3.6.
- Moved scripts to sub-module and updated entry points.
- Switched to Poetry for dependency management.
- Using black for format consistency.

### Fixed
- Analysis of 2D images now works properly.
- Dependency conflicts.



## [3.3.5] - 2019-02-13
### Added
- `czi_to_tiff`
    + Support for 2D conversion.
- `nd2_to_tiff`
    + Support for 2D conversion.
- `tiff_findoof`
    + Parallelization.

### Changed
- Now using `ggc` functions to export settings in: `gpseq_anim`, `gpseq_fromfish` and `tiff_auto3dseg`.
- `czi_to_tiff`
    + Refactored for easier development.
- Now using `seaborn` for color palettes.

### Fixed
- Missing dependencies in `setup.py`.
    + Installation through pypi.org.
- `tiff_auto3dseg`
    + `--neighbour` now works properly.
    + Provided better error message when 2D mask folder not found.
- `nd2_to_tiff`: case of single-channel stacks.
- `gpseq_anim`: fixed skipping of boxplot step (previously broke skipping)



## [3.3.4] - 2018-09-21
### Added
- `czi_to_tiff`: to convert CZI files to TIFF.
- `gpseq_fromfish_merge`
    + Option for no-date prefix to output.
    + Support for date in dataset name.
- `tiff_findoof`
    + Silent mode with `-s`.

### Changed
- `gpseq_anim`
    + Now using `ggc.check_threads()`.
- `tiff_auto3dseg`
    + Now using `ggc.check_threads()`.
- Fixed `--compressed` option label (now compatible with ImageJ).
- `tools.plot.save_tif`
    + Added support to retain voxel resolution in TIFF metadata (ImageJ compatible).
- `nd2_to_tiff`
    + Now saves Resolution (XYZ) metadata when exporting.



## [3.3.3] - 2018-09-11
### Added
- `tiff_split v1.1.0`
    + Allowed for splitting into non-square/cubic images.
    + Option to change splitting orientation.



## [3.3.2] - 2018-08-28
### Fixed
- `gpseq_fromfish`
    + Incompatibility with newer version of `tiff_auto3dseg` (`v3.1.0+`) script caused by unexpected channel axis in input mask.

### Added
- `gpseq_fromfish`
    + Enforcing re-slicing to 3 axes of the input masks, to match the input images.
    + Readable error message in case of inconsistent shape between input mask and image, reverting to binarization in that case.



## [3.3.1] - 2018-08-23
### Fixed
- `gpseq_fromfish`
    + Bug in nuclear semi-axes length calculation.



## [3.3.0] - 2018-08-22
### Added
- Clearer documentation for homologue copy pairs to `gpseq_fromfish` and `gpseq_fromfish_merge`.
- `gpseq_fromfish`
    + Additional help page with `-H`.
    + Option for 0-indexed input.
    + New columns to output nuclear table:
        * Using `slice`, `row` and `col` for coordinates, instead of `z`, `y` and `x`.
        * `box_start_slice`, `box_start_row`, `box_start_col`: nuclear box starting point 1-indexed coordinates (integer).
        * `box_end_slice`, `box_end_row`, `box_end_col`: nuclear box ending point 1-indexed coordinates (integer).
        * `com_slice`, `com_row`, `com_col`: nuclear mask center of mass 0-indexed coordinates (float).
    + New columns to output compartment table:
        * `*_slice_component`, `*_row_component`, `*_col_component`: with the components of the three major nuclear axes over the image dimensions.

### Changed
- Clarified warning when input image axes do not match with metadata.
- `gpseq_fromfish`
    + Split script help page in `-h` for attributes and standard help, and `-H` for more descriptive and readable text.
    + FISH coordinates can now be floating point (integer is not enforced anymore). Lamina/Center distances are interpolated on the regulard grid.
    + Silenced low contrast warnings when saving stacks in debugging mode.
    + Changed extension of output tables to `.tsv`, for consistency with actual formatting.
    + Output compartment table renamed to `nuclear_compartment.out.dilate*.*`.
- `gpseq_fromfish_merge`
    + Clarified `no copy pairs found` warning message.
    + `--aspect` default changed to `300. 130. 130.` (ZYX).
    + Now compatible with both new and old nuclear compartment table naming.
- `tiff_findoof`
    + Gradient magnitude mode now default, switch to intensity sum mode with `-S` or `--intensity-sum`.

### Fixed
- `gpseq_anim`
    + Crash when 2D mask folder is not provided.
    + Wrong variable name in `anim.series.Series`.

### Removed
- `gpseq_fromfish`
    + Removed `com` (center of mass) column from the output dot table. Now the same information is available in the nuclei table (although the CoM coordinates are box-wise, and not image-wise).
    + Exporting dot table before `Allele` (homologue copy pair) labeling.



## [3.2.1] - 2018-08-20
### Fixed
- Now writing `tiff` files with proper `axes` metadata.

### Changed
- `gpseq_fromfish`
    + Removed double-negation in plot settings confirmation.



## [3.2.0] - 2018-08-16
### Added
- `tiff_auto3dseg`
    + Added option to discard objects touching the Z borders during segmentation.
    + `-F` option for dilate-fill-erode operation (10 as default).
- `tiff_findoof`
    + `-G` option for gradient magnitude mode.
    + `-R` rename option.

### Fixed
- `gpseq_fromfish`
    + Allowed for missing labels.
    + Removing `NaN` and `Inf` when plotting aggregated visualization.
    + Bug due to high imprecision when calculating angles.
- `tiff_auto3dseg`
    + Combination with 2D mask now moved before border clearing.



## [3.1.0] - 2018-07-16
### Added
- `tiff_auto3dseg`
    + Added option to combine 2D mask with 3D ones.
    + Allowed for labeled input (2D masks).
- `gpseq_fromfish` & `gpseq_anim`
    + Added option for 2D-to-3D mask combination (also labeled input).

### Changed
- `tiff_auto3dseg`
    + Removed superfluous `nargs` argument in `add_argument`.



## [3.0.4] - 2018-07-02
### Added
- Version info to package dependencies in `setup.py`.
- More flexible adaptive threshold behavior in `tools.binarize`.
- Option to use dilated mask only for dots assignment, and not for distance calculation, as `gpseq_fromfish --dots-assignment-only`.
- (Un)Compression log to `tiffcu`.

### Changed
- Now depending on `scikit-image v0.14.0` (upgraded package).
- Output directory now optional in `nd2_to_tiff`.
- `gpseq_fromfish` now not crashing when a tiff file is corrupted, only an error message is showed.
- Better error message during initialization if a series does not have either DNA or Signal channel(s).

### Fixed
- Minor bugs in `binarizer` that crashed the script with certain parameter combinations.
- Local threshold behavior.
- Order of output table columns in `gpseq_fromfish`.
- `nd2_to_tiff` now compatible with single-channel nd2 files.
- `nd2_to_tiff` now exporting tiff files with proper axes order.
- Mask import in `gpseq_anim`.



## [3.0.3] - 2018-06-29
### Fixed
- Typo in `gpseq_fromfish` that rendered the `--dilate-Z` flag useless.



## [3.0.2] - 2018-06-22
### Fixed
- Bug that crashed `gpseq_anim` when option `-m` is not used.



## [3.0.1] - 2018-06-21
### Fixed
- `tiff_auto3dseg` minor bug due to unbound `staticmethod`.



##  [3.0.0] - 2018-06-12
### Added
- Enzyme diffusion simulation.
- `gpseq_fromfish`
    + Settings export.
    + Debug mode.
- `gpseq_anim`
    + `--do-all/-y` option.
    + `-M`, `-m`, `--labeled` and `--compressed` to export/import masks.

### Changed
- Moved distance functions to separate `tools` sub-module.
- Now using `ggc` to confirm settings in `gpseq_anim`.

### Fixed
- `return` bug crashing the pipeline when a condition has no nuclei.



## [2.1.8] - 2018-06-07
### Fixed
- `datetime` import in `gpseq_anim`.

### Added
- Volume density merge through `gpseq_fromfish_merge`.



## [2.1.7] - 2018-06-01
### Added
- Nuclear density merge through `gpseq_fromfish_merge`.
- Volume profile to `gpseq_fromfish`.

### Fixed
- Center as percentile option is now used also for the nuclear density calculation in `gpseq_fromfish`.



## [2.1.6] - 2018-05-31
### Fixed
- Bug in `gpseq_fromfish` due to wrong column labeling.
- Bug in `gpseq_fromfish_merge`.

### Added
- Python version, command line and timestamp to `gpseq_anim` exported settings file.



## [2.1.5] - 2018-05-07
### Added
- Option `--min-radius`/`-R` to `gpseq_anim` to change the XY size filter of segmented objects.



## [2.1.4] - 2018-05-03
### Added
- Intensity sum on sum Z-projection for each nucleus as `flat_sumI`.
- Density profile for each nucleus to `gpseq_anim` and `gpseq_fromfish`.
- Settings confirmation to `gpseq_fromfish`



## [2.1.3] - 2018-04-24
### Added
- Plotting option to `tiff_findoof` script.
- Option to skip hole filling in segmented masks (`--no-hole-filling`) in `gpseq_anim`.
- Option to define nucleus center as top percentile (`--center-as-percentile`).
- Settings to file in `gpseq_anim`.

### Fixed
- Axes order in `nd2_to_tiff`.
- Lamin distance calculation in 3D.



## [2.1.2] - 2018-04-10
### Fixed
- Fixed dataset aggregated view file names.
- Extension now present in `nd2_to_tiff` output.



## [2.1.1] - 2018-03-28
### Changed
- Z dilation is now turned off by default.

### Fixed
- Compartment assignment.
- Fixed colors in aggregated visualization.



## [2.1.0] - 2018-03-28
### Added
- Aggregated visualization per field of view and for the whole dataset, per channel and for all channels.
- Dilation can be enforced to XY dimensions only (`--no-Z-dilation`).

### Changed
- In 3-ortho-view plot nuclei are now occupying as much space as possible.
- In 3-ortho-view plot are now colored based on the channel.

### Fixed
- Now rotating only on XY plane for nuclear compartmentalization.

### Removed
- Option `--annotate-compart` is now the default behaviour.



## [2.0.1]
### Fixed
- Fixed sub packages import.
- Now possible to install not in development status.



## [2.0.1]
### Added
- `tiff_split` from `tiff-tools-gg`.
- `--skip-channels` option to `gpseq_fromfish`.
- `nd2_to_tiff` conversion script.
- Progress bar to `tiff_findoof`.

### Fixed
- Fixed `gpseq_fromfish` missing `istruct` passing.
- Input directory list issue in `gpseq_fromfish_merge`.
- Allowed for only partial information in `gpseq_fromfish_merge` input.
- TIFF mask import/export in `gpseq_fromfish`.



## [2.0.0]
### Changed
- Cleaned up `const.py`.
- `wraps` module renamed `anim`, will contain all wrappers for GPSeq standard image dataset analysis.
- `main` library moved to `anim` module.
- Cleaned up and re-structure comment style in all libraries and modules.

### Added
- License text in `bin/gpseq_anim`.
- TIFF automatic 3D segmentation script from `tiff-tools-gg`.
- TIFF (un)compress script from `tiff-tools-gg`.
- Out Of Focus detection/plot scripts from `tiff-tools-gg`.
- FISH lamin distance calculation script from `dotter2gpseq`.
- `gpseq_fromfish` output merge from `dotter2gpseq`.



## [1.1.0] - 2018-03-13
### Added
- Capturing TIFF read issues.
- Voxel aspect unit of measure.
- `Max` profile.

### Fixed
- Now borders are properly cleared when analyzing in 3D.
- Nuclear threshold plot size.
- Report format.
- Plot titles.

### Changed
- Cleaned package structure.
- Renamed `pygpseq-run.py` to `gpseq_anim` and allowed for easy installation.
- `VERSION` constant now source only from `setup.py`.
- Setup base for unit testing with `nosetests`.
- Sigma used for smoothing and density curve production are now separate.



## [1.0.0]
### Changed
- Now normalizing lamin distance over its sum with central distance.

### Added
- Sigma option to `pygpseq-run.py`.


## [0.1.2]
### Changed
* CSV output files in separate folder.
* X-axis intercepts of profiles and profiles 1st/2nd derivatives are reported.
* Profile areas are reported.
* Single-nucleus data summaries are exported in CSV format.
* Background levels are now reported as barplots.

### Added
* Automatic stop if no conditions are found.
* Automatic stop if no series are found in a condition.



## [0.1.1]
### Changed
- Better documentation is now available.
- DocStrings are now in Napoleon format.
- Improved report.
- Nuclear selection can now be skipped.
- Nuclear selection can now be perform with up to 3 features.
- Moved image binarization in a separate class.
- Generalized channel reading function.

### Added
- Distribution comparison using Wilcoxon-Mann-Whitney U test.
- Implemented different middle-section definitions.



## [0.1.0] - 2017-07-26




[Unreleased] https://github.com/ggirelli/gpseq-img-py
[3.4.1] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.4.0
[3.3.5] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.5
[3.3.4] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.4
[3.3.3] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.3
[3.3.2] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.2
[3.3.1] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.1
[3.3.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.3.0
[3.2.1] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.2.1
[3.2.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.2.0
[3.1.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.1.0
[3.0.4] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.0.4
[3.0.3] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.0.3
[3.0.2] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.0.2
[3.0.1] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.0.1
[3.0.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v3.0.0
[2.1.8] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.8
[2.1.7] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.7
[2.1.6] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.6
[2.1.5] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.5
[2.1.4] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.4
[2.1.3] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.3
[2.1.2] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.2
[2.1.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.1.0
[2.0.1] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.0.1
[2.0.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v2.0.0
[1.1.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v1.1.0
[1.0.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v1.0.0
[0.1.3] https://github.com/ggirelli/gpseq-img-py/
[0.1.2] https://github.com/ggirelli/gpseq-img-py/
[0.1.1] https://github.com/ggirelli/gpseq-img-py/
[0.1.0] https://github.com/ggirelli/gpseq-img-py/
