# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased
### Added
- Capturing TIFF read issues.
- Voxel aspect unit of measure.

### Fixed
- Now borders are properly cleared when analyzing in 3D.
- Nuclear threshold plot size.
- Report format.

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
[1.0.0] https://github.com/ggirelli/gpseq-img-py/releases/tag/v1.0.0
[0.1.3] https://github.com/ggirelli/gpseq-img-py/
[0.1.2] https://github.com/ggirelli/gpseq-img-py/
[0.1.1] https://github.com/ggirelli/gpseq-img-py/
[0.1.0] https://github.com/ggirelli/gpseq-img-py/
