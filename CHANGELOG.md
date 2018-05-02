# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased
### Added
- Intensity sum on sum Z-projection for each nucleus as `flat_sumI`.



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
