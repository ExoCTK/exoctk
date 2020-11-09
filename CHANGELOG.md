# Changelog

This file tracks all major changes in each `exoctk` release.

## [1.1.1] - Unreleased

### Added

- Added `CHANGELOG.md` (this file!) to repo

### Changed

- Updated readthedocs/JDox links
- `pyup` updates to several packages

## [1.1.0] - 2020-10-12

### Added

- Added Jupyter notebook for Phase Constraint tool
- Added ACES and ATLAS9 model grids as child classes for Limb Darkening tool

### Fixed

- Fixed MIRI trace shape in Contamination and Visibility tool
- Fixed table truncation in Limb Darkening Calculator results table

## [1.0.0] - 2019-12-31

### Added

- Added Amazon AWS cloud computing support
- Added Contamination and Visibility Tool support for all JWST instruments
- Added Phase Constraint Tool

## [0.4.0] - 2019-08-09

### Added

- Added tests for atmospheric retrieval code
- Added form validation to Web app with `flask_wtf`
- Added Lightcurve Fitting tool
- Added generic model grid to forward models
- Added `jwst_gtvt` package dependency for Contamination and Visibility tool
- Added planetary parameter support with `exomast`
- Added `platon` support to forward models
- Added dependency support with `pyup`

### Removed

- Package data in favor of package data download

## [0.2.0] - 2018-11-11

### Added

- PEP257 and PEP8 compliance for all tools
- Updated LDC plotting and notebook
- Documentation support with `sphinx-build`

### Changed

- Changed `ldc` function to `LDC` class
- Renamed top level directory `exoctk` from `ExoCTK`
- Major `README` updates

### Removed

- Built-in filter support removed in favor of `svo_filters` package dependency
- Removed `astropy_helpers` support since it was not helping

## [0.1.4] - 2017-07-10

### Added

- Added Contmination and Visibility tool
- Added additional filter support to LDC including NIRISS GR700XD and custom filter profile

## [0.1.3] - 2017-07-07

### Added

- Line Forward Models support with `chimera`

## [0.1.2] - 2017-06-30

### Added

- Added precomputed pickle file of bandpass throughputs for speedier LDC calculations

## [0.1.1] - 2017-06-29

### Fixed

- Minor bug fixes to initial release

## [0.1.0] - 2017-06-28

### Added

- Limb darkening calculator (LDC)
- Forward modeling with exotransmit (PAL)
- Transit observation tools using PandExo (TOT)
- Groups and integrations calculator (TOR)
- Model grid handling
- References tracking
- Jupyter notebooks for LDC, PAL, and TOT modules
- SVO Filter Profile Service for photometric bandpasses

[1.1.1]: https://github.com/ExoCTK/exoctk/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/ExoCTK/exoctk/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/ExoCTK/exoctk/compare/v0.4.0...v1.0.0
[0.4.0]: https://github.com/ExoCTK/exoctk/compare/v0.2.0...v0.4.0
[0.2.0]: https://github.com/ExoCTK/exoctk/compare/v0.1.4...v0.2.0
[0.1.4]: https://github.com/ExoCTK/exoctk/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/ExoCTK/exoctk/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/ExoCTK/exoctk/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/ExoCTK/exoctk/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/ExoCTK/exoctk/releases/tag/v0.1.0