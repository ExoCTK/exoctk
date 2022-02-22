# Changelog

This file tracks all major changes in each `exoctk` release.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.3] - 2022-02-01

### Fixed

- Date format in hover tool of contam visibility plots
- Logging issue with production server.

### Changed

- Removed support for "uniform" (constant) limb-darkening law.
- Removed lmfit from requirements.

## [1.2.2] - 2021-07-09

### Added

- Sweep to update code to match PEP8 standards.
- Extra authors on citation information to match current working DOI.

### Fixed

- Phase-constraint bug that didn't change `eccentricity` to `nan` when not found by the target resolver.

## [1.2.1] - 2021-06-09

### Added

- NIRSpec grism support for the limb darkening calculator
- Test module for the `forward_models.py` module

### Fixed

- Broken forward modeling page of the Web application

## [1.2.0] - 2021-05-10

### Added

- Citation information in the `README` file (see PR #472).
- Contribution guide on the project's wiki (https://github.com/ExoCTK/exoctk/wiki/Git-Workflow).
- Procedures for ExoCTK release (https://github.com/ExoCTK/exoctk/wiki/Release-Procedure).
- Support for installation via PyPI (see PR #385).
- JWST instrument throughputs for limb-darkening calculations (see PR #458).

### Fixed

- Dropped support for `python3.6` (see PR #471).
- Installation problem with Python 3.8 (see issue #465 and PR #469).
- Pip & asteval enviornment issues (see PR #468).

### Changed

- Missing links to JDox on contamination tool (see PR #433).
- ExoCTK data is now downloaded as-needed for individual tools (see PR #488).

## [1.1.1] - 2021-05-10

### Added

- Added `CHANGELOG.md` (this file!) to repo

### Fixed

- WFC3_IR.G102 filter in LDC

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

[1.2.1]: https://github.com/ExoCTK/exoctk/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/ExoCTK/exoctk/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/ExoCTK/exoctk/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/ExoCTK/exoctk/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/ExoCTK/exoctk/compare/v0.4.0...v1.0.0
[0.4.0]: https://github.com/ExoCTK/exoctk/compare/v0.2.0...v0.4.0
[0.2.0]: https://github.com/ExoCTK/exoctk/compare/v0.1.4...v0.2.0
[0.1.4]: https://github.com/ExoCTK/exoctk/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/ExoCTK/exoctk/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/ExoCTK/exoctk/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/ExoCTK/exoctk/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/ExoCTK/exoctk/releases/tag/v0.1.0
