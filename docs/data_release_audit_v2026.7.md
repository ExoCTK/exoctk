# ExoCTK v2026.7 external-data release audit

## Executive summary

This is an audit and release plan for the full ExoCTK v2026.7 release at
`origin/main` commit `d5bd25ed5d293abbe21a8e38aa490a1acb516219`. It is not a
Pandeia-only data refresh, and it does not implement a production migration.

The current main branch still exposes the Fortney and generic forward models,
groups/integrations, ACES and ATLAS9 limb darkening, and contamination and
visibility tools. The old data cannot therefore be retired as a unit. The
evidence supports independently versioned components, not another monolithic
archive.

The principal findings are:

* `groups_integrations.json` is no longer the canonical runtime table. Current
  code reads the packaged
  `exoctk/data/groups_integrations/groups_integrations_input_data.json`.
  That small table must nevertheless be regenerated/revalidated because its
  Pandeia provenance and builder are missing and current main added NIRSpec
  PRISM and NIRCam readout behavior.
* The Fortney (46.5 MB) and generic (2.31 GB) databases remain consumed by
  public website routes and Python functions. Their on-disk schemas still
  match current loaders. Keep them unchanged as optional forward-model data;
  do not reconstruct their uncertain scientific builders for this release.
* ACES (16.82 GB locally) and ATLAS9 (64.68 MB locally) remain required for the
  supported limb-darkening feature. Pandeia's Phoenix/Kurucz SED APIs return
  disk-integrated wavelength/flux spectra, not the wavelength-by-`mu`
  intensities that `ModelGrid` and `LDC` require. They are not replacements.
* Current contamination calculations require SOSS trace FITS, NIRISS order-0
  NumPy templates, and newly exposed NIRCam DHS F322W2/F444W parent trace
  templates. The directly inspected Box contamination ZIP has SOSS and an old
  DHS-F150W2 family, but no F322W2/F444W directories. Candidate DHS files exist
  in the local old-data tree but have no established publication provenance.
* No contamination product should be regenerated or released until PR #725
  and PR #726 are merged (or their final behavior is otherwise on main) and
  order-0 calibration receives explicit scientific approval. Current code
  multiplies order 0 by `1.5e3`; the absolute calibration has not been
  validated. This audit neither changes source classification nor Gaia TAP.
* Target-specific contamination HDF5 databases and `exoctk_log.db` are mutable
  runtime state. They should be initialized/migrated and invalidated by
  provenance, not shipped as release data.
* The 170+ packaged JWST throughput tables trace to Pandeia refdata v1.5.2.
  Current Pandeia exposes `InstrumentFactory.get_total_eff`, but a complete
  configuration, units, wavelength, and numerical equivalence test has not
  been performed. Regenerate or replace them only in a focused follow-up.
* A temporary build from an unmodified source export produced a 181,432-byte
  wheel with 54 files and a 168,964-byte sdist. Neither contained any
  `exoctk/data`, Flask template, or static file. `setup_package.py` defines a
  package-data list but `setup.py` does not use it. This is a release blocker.
* `utils.py` constructs eight `...v1.2.tar.gz` STScI URLs. All eight return
  `file_not_found`. The corresponding unversioned STScI URLs used in the old
  data Dockerfile are accessible. Website CI instead downloads eight Box URLs,
  including one duplicate contamination URL and an `all.zip` that contains
  only an empty directory.
* The clean scientific portion of the local old installation is approximately
  23.28 GB. A leftover 10.52 GB `output.zip` raises actual disk use to roughly
  33.8 GB. A future component total cannot be stated until contamination and
  groups products are regenerated, but optional separation avoids imposing
  the 16.82 GB ACES grid on all users.

The machine-readable conclusions are in
`docs/data_release_inventory_v2026.7.json`.

## Decision table

| Data family | Current source | Approximate size | Current consumers | Feature status | Provenance | Generator | Dependency replacement | v2026.7 decision | Required follow-up |
|---|---|---:|---|---|---|---|---|---|---|
| Packaged groups table | `exoctk/data/groups_integrations/...json` | 253 KB | Groups API and `/groups_integrations` | Active | Pandeia version unknown | Not found | Direct Pandeia equivalence unresolved | Regenerate | Recover/write versioned builder |
| External groups table | Box/STScI/local `groups_integrations.json` | 523 KB | Stale tutorial only | Superseded | 2020, engine unknown | Not found | Packaged table is canonical | Remove | Docs/package cleanup |
| Fortney DB | Box/STScI/local | 46.50 MB | `fortney_grid`, `/fortney` | Active | Fortney family; builder unknown | Not found | None | Keep | Optional manifest/provenance |
| Generic DB | Box/STScI/local | 2.313 GB | `generic_grid`, `/generic` | Active | Unknown | Not found | None | Keep | Optional manifest/provenance |
| ACES intensities | STScI/local | 16.824 GB | ModelGrid/LDC/website | Active | PHOENIX ACES filenames/headers | Upstream unknown | Pandeia SED is not `mu`-dependent | Keep | Separate limb-darkening component |
| ATLAS9 intensities | STScI/local | 64.68 MB | ModelGrid/LDC/website | Active/default | Converted ATLAS9 | `helpers.convert_ATLAS9` | Pandeia SED is not `mu`-dependent | Keep | Record original inputs |
| ModelGrid indexes | Beside raw grids | Unknown complete size | ModelGrid cache loader | Active accelerator | Path-dependent/partial locally | `ModelGrid` | None | Regenerate | Versioned, complete, path-neutral cache |
| SOSS traces | Box/STScI/local | 1.260 GB local SUBSTRIP256 | Contamination worker/API | Active | Pandeia product, version absent | `make_contam_traces.py` | Pandeia is generation input only | Regenerate | Wait for blockers; add manifest |
| NIRCam DHS traces | Local candidate only | 63.60 MB | F322W2/F444W contamination | Newly active | Unknown | Partial generator | Pandeia is generation input only | New | Validate/regenerate after blockers |
| NIRISS order 0 | Box/local | 61 KB clean local | SOSS contamination | Active | Unknown; operational scaling | Not found | None | Regenerate | Scientific approval first |
| Other legacy traces | Local mixed tree | Unknown clean subset | Broad Python aperture path; not web contamination | Ambiguous | Mixed | Partial | Pandeia input only | Conditional regenerate | Decide API support first |
| Contamination HDF5 | `$EXOCTK_CONTAM_CACHE` | Deployment-dependent | Runtime cache | Active/mutable | Live queries/templates; not recorded | `precompute.py` | N/A | Runtime | Version/invalidate/migrate |
| Logging DB | `$EXOCTK_DATA/exoctk_log` | 28–90 KB observed | Website logging | Active/mutable | Schema in code | `create_db` | N/A | Runtime | Stop distributing populated DB |
| JWST throughputs | Python package | 4.46 MB | Throughput/LDC | Active | Pandeia refdata v1.5.2 | `generate_JWST_throughputs` | Candidate, not proven | Regenerate | File-by-file 2026.7 equivalence |
| Small contamination tables | Python package | 217 KB family | Gaia temperature, wavelength, visibility, precompute | Active/mixed | Partly unresolved | Mixed | `jwst_gtvt` replaces only newer visibility path | Keep | Document/split runtime vs generation |
| Unused 2MASS/long ephemeris | Python package | Small | None found | Unexposed | Legacy | Unknown | N/A | Remove | Focused cleanup after final search |
| Packaged test grids/cache | Python package/tests | 6.52 MB | Unit tests | Fixtures | Partial | Existing code | N/A | Keep | Label and add provenance |
| macOS metadata | Published ZIPs/local | Archive-dependent | None | Artifact | Finder/Archive Utility | N/A | N/A | Remove | Exclude from new archives |

## 1. Scope and methodology

The audit branch is based directly on the latest fetched `origin/main`. No
contamination PR was merged or copied. The audit used:

* repository-wide exact-name, environment-variable, URL, route, loader, test,
  documentation, and history searches;
* a comparison from the latest public release tag to `origin/main`;
* metadata-only HTTP probes that follow redirects with a one-byte range GET
  fallback;
* direct read-only inspection of downloaded archives in a temporary directory;
* read-only size, FITS, HDF5, JSON, pickle-location, and SQLite schema checks on
  `/Users/tbell/exoctk_data`;
* a temporary wheel/sdist build from `git archive HEAD`; and
* source inspection of installed Pandeia 2026.2 with local 2026.7 JWST
  refdata. This mixed installed state is evidence about API/data shape, not a
  claim that Pandeia engine 2026.2 is the release engine.

No scientific product was generated, no archive was extracted into the old
data tree, and no file under `/Users/tbell/exoctk_data` was changed.

## 2. Baselines examined

### Latest previous public release

GitHub's latest published release is `v1.2.6.3a`, tag commit
`2aee86159ea06b736d24d0cbb72b242e68dfd902`, published 2025-09-25. It is an
ancestor of current main.

### Current main

The audit base is `d5bd25ed5d293abbe21a8e38aa490a1acb516219` (2026-07-20),
543 commits after the tag by `git describe`. Material post-release work
includes precomputed/living contamination caches, NIRCam DHS contamination,
NIRSpec PRISM and NIRCam readout handling in groups/integrations, local
Pandeia-RC Docker support, order-0 templates/scaling, and expanded website
tests.

### Deployment configuration

`docker/compose.yaml` bind-mounts `$EXOCTK_DATA`, a separate
`$EXOCTK_CONTAM_CACHE`, and read-only Pandeia refdata. ExoCTK and worker images
receive `$EXO_VERSION`, but no production value is tracked. The Pandeia Docker
service currently declares 2026.2. GitHub's deployments API returned an empty
list; this repository contains build/test workflows but no production
deployment workflow.

### Public website

On 2026-07-20 local time, `https://exoctk.stsci.edu/` returned HTTP 200 and a
footer link labeled `exoctk v1.2.6`. The same label is hard-coded in current
main (`templates/base.html:149`), so it is not evidence of an exact deployed
commit. The page exposes no Git SHA, data manifest, or cache version.

### Published archives and local installation

Both Box and STScI URL families were mapped separately. The local installation
was treated as an independent older baseline, not as proof of publication or
current need.

## 3. Uncertainty about the deployed website version

No direct evidence ties together an exact deployed Git SHA, Python package
version, external-data archive set, Pandeia/refdata version, and contamination
cache. The public footer supports only the statement that the site presents
itself as v1.2.6. Current main retains that same hard-coded footer. GitHub has
no deployment records for the repository, and production `.env` values are not
tracked. This uncertainty must remain in release notes until deployment emits
a build/data/cache manifest or an operator records the live container state.

## 4. Current external-data consumers

### Import-time behavior

`utils.py:62` chooses `$EXOCTK_DATA` or `~/exoctk_data`. Outside CI/RTD it
creates missing family directories at import time (`utils.py:74-90`). That is
an existence/setup side effect, not proof each family is consumed.

`app_exoctk.py:63-70` opens or creates `exoctk_log/exoctk_log.db` while the web
module imports. `field_simulator.py:175` loads the packaged Gaia color table at
import. `contamination_figure.py:31` resolves the packaged wavelength map at
import.

### Runtime opens

* Fortney: `forward_models.py:67-72` checks the family and opens
  `fortney_models.db` with SQLAlchemy.
* Generic: `forward_models.py:210-216` checks the family and opens
  `generic_grid_db.hdf5`; `rescale_generic_grid` reads catalog, wavelength, and
  the selected spectrum.
* Groups/integrations: four functions open the caller's JSON. The website and
  tests pass the packaged input table, not `$GROUPS_INTEGRATIONS_DIR`.
* Limb darkening: `ACES` and `ATLAS9` build paths under `$EXOCTK_DATA`.
  `ModelGrid` opens raw FITS intensity, `mu`, abundance, and wavelength arrays
  and optional pickle/HDF5 accelerators.
* Contamination: `field_simulation` checks `exoctk_contam`; `get_trace` globs
  mode/temperature FITS; `get_order0` loads NIRISS `.npy` templates. The
  optional HDF5 cache is opened in append mode and may be created on a miss.
* Logging: the web application writes form submissions to the SQLite DB.

No current CLI entry points are declared in `pyproject.toml`. Website workers
consume groups and contamination data through Celery tasks.

## 5. Box archive mapping

Website CI references eight Box download commands but only seven unique URLs;
`zhrx...` (contamination) appears twice.

| URL ID | Resolved file | Compressed | Members / expanded | SHA-256 | Inspection and classification |
|---|---|---:|---:|---|---|
| `grztwf220jnvamzv8ugb7phj9tlwnoi1` | `all.zip` | 162 B | 1 empty directory / 0 B | `fc8b04...1045` | Direct; obsolete placeholder |
| `zhrx2r6kgeovcnj78d4wicb7m2009qca` | `exoctk_contam.zip` | 2,837,723,039 B | 648 / 3,797,697,305 B | `e0e13a...1045` | Direct; partially required and stale/incomplete |
| `4q75kt2xelk3zox26vii1i3clnm7nzt7` | `exoctk_log.zip` | 8,581 B | 4 / 90,536 B | `e5bce1...3a99` | Direct; mutable runtime state |
| `702r4nkccc484noiggq5xt13agbkkevp` | `fortney.zip` | 26,998,140 B | 4 / 46,498,216 B | `536f90...776` | Direct; required unchanged |
| `p64c53rky00ojz1gvabe9hqirpat9azz` | `generic.zip` | 1,815,267,582 B | 4 / 2,313,342,496 B | `2ea194...9e9e` | Direct; required unchanged |
| `1ag5x8ii1ko8h6lnoxwxljgdfw2mp894` | `groups_integrations.zip` | 112,286 B | 6 / 524,220 B | `2d2ef5...acef` | Direct; external copy superseded |
| `vcn1ff9nzg2qx3pnkmsz17dujjdzlsh1` | `modelgrid.zip` | 10,518,396,973 B | Remote central directory not fetched | Unknown | Metadata-only; local `output.zip` has identical byte size but identity is unresolved |

The directly inspected contamination ZIP contains only these scientific trace
roots: `NIRISS` (with macOS sidecars), `NIS_SUBSTRIP256`, `NIS_SUBSTRIP96`, and
`NRCA5_DHS_F150W2`, plus order-0 files. It does not contain the current
F322W2/F444W parent directories. It has 327 macOS metadata members. The generic
and Fortney ZIPs each have two metadata members. The groups ZIP has four. No
inspected ZIP has duplicate or unsafe member paths.

The local `output.zip` has the same 10,518,396,973-byte size advertised for
Box `modelgrid.zip`, 2,785 members, 17,218,505,022 expanded member bytes, 1,393
macOS metadata members, and SHA-256
`73355d9645cdb2e21bb0edd30846432e6673a3a6ec555fe29f84f6d9cd470d82`.
Without a remote checksum this is not claimed to be the Box object.

## 6. STScI archive mapping

`utils.py` derives `PATCHVER = v1.2` from `_version.py` 1.2.6.4 and constructs
eight versioned STScI URLs. Every versioned URL returned a 404 JSON
`file_not_found` response. Therefore `download_exoctk_data()` is currently
broken independently of any data-content decision.

The old Dockerfile references accessible unversioned archives:

| Archive | Compressed size | Direct inspection | Contents / decision |
|---|---:|---|---|
| `exoctk_contam.tar.gz` | 1,350,834,690 B | Metadata only | Required but stale/incomplete; direct Box counterpart inspected |
| `groups_integrations.tar.gz` | 110,821 B | Direct, SHA `389b89...b0af` | One 523,408-B JSON; same payload SHA as Box/local; superseded |
| `fortney.tar.gz` | 26,997,280 B | Direct, SHA `8ba228...0e7d` | One 46,497,792-B DB; same payload SHA as Box/local; keep |
| `generic.tar.gz` | 1,816,705,192 B | Metadata only | Required unchanged; Box counterpart directly inspected |
| `exoctk_log.tar.gz` | 862 B | Direct, SHA `c412a0...0d4` | 28,672-B DB; mutable state; remove from release data |
| `modelgrid_ATLAS9.tar.gz` | 70,913,836 B | Direct, SHA `dcdca6...db2c` | 333 members, 87,167,542 expanded bytes, no macOS/unsafe paths; keep optional |
| `modelgrid_ACES_1.tar.gz` | 5,263,356,499 B | Metadata only | Required optional ACES part 1 |
| `modelgrid_ACES_2.tar.gz` | 5,260,200,291 B | Metadata only | Required optional ACES part 2 |

The STScI tar payload hashes for Fortney
(`24be9157...d504`) and groups (`a6cde638...896d`) match their Box and local
counterparts. The STScI log DB differs from Box and local log DBs, confirming
that it is mutable/version-divergent state.

## 7. Local old-data inventory

The read-only tree contains:

| Family | Clean file bytes | Notable content |
|---|---:|---|
| `exoctk_contam` | 4,033,067,841 | 21 trace directories, order 0, legacy IDL/Python artifacts |
| `modelgrid` | 16,888,679,683 | 1,056 ACES scientific FITS, 331-ish ATLAS/index files, derived caches |
| `generic` | 2,313,342,072 | One HDF5 database |
| `fortney` | 46,497,792 | One SQLite database |
| `groups_integrations` | 523,408 | Canonical old JSON plus AppleDouble sidecar |
| `exoctk_log` | 90,112 | Mutable SQLite DB |
| Clean family total | 23,282,200,908 | About 23.28 GB decimal |
| Leftover `output.zip` | 10,518,396,973 | Modelgrid archive candidate, not installation data |

The tree also has `__MACOSX`, `.DS_Store`, `._*`, two tiny JSON error bodies
misnamed `.tar.gz`, and the large leftover ZIP. Presence alone was not used to
mark any family required.

The local contamination tree is broader than the Box contamination ZIP:
MIRI LRS, old NIRCam names, NIRCam grism filters, nine NIRSpec BOTS modes,
NIRISS legacy products, SOSS subarrays, and DHS-F150W2. This mixed tree cannot
be assigned a single provenance or support status.

## 8. Files packaged directly with ExoCTK

`exoctk/data` has 211 tracked files and 23,976,915 bytes:

* current groups/integrations JSON;
* Gaia color, SOSS wavelength, target-list, ephemeris, and legacy 2MASS tables;
* 170+ JWST/TESS/CHEOPS throughput tables;
* eight ATLAS9 test FITS plus test indexes;
* FITS template, bibliography, Vega table, images, and directory sentinels.

The package-data build check is negative: current wheel and sdist contain none
of these files and no Flask templates/static files. `include-package-data =
true` is insufficient without selected data in the sdist, and the
`get_package_data()` function is unused. Source-checkout tests can hide this
because imports resolve the checkout's `exoctk` directory.

The `ModelGrid_tmp.fits` template and bibliography are generation/reference
inputs. The small model grid and contamination HDF5 are fixtures. Images are
documentation/UI assets, not external scientific data.

## 9. Data replaceable by dependencies

No external scientific family is proven safely replaceable in this audit.

Pandeia 2026.2's `pandeia.engine.sed.Phoenix` calls
`stsynphot.grid_to_spec('phoenix', teff, metallicity, log_g)` and returns
wavelength in microns and flux in mJy. The local 2026.7 refdata SED config
describes those disk-integrated spectra. ACES and ATLAS9 limb-darkening readers
instead require intensity shaped by atmospheric parameters, angle/`mu`, and
wavelength. The sample ACES primary array is `[78, 25500]` and ATLAS9 is
`[17, 1217]`. Scientific meaning and dimensionality are different.

Pandeia does replace copied stellar SEDs as an input to contamination trace
generation: `make_contam_traces.py` already uses `sed_type='phoenix'`,
`build_default_calc`, and `perform_calculation`. It does not replace the
generated detector trace FITS that current runtime code opens.

The throughput files are a plausible installed-dependency replacement because
`InstrumentFactory.get_total_eff` is the existing generator API and produces
micron plus dimensionless PCE. Replacement is not approved until every file's
configuration mapping, coverage, sampling, units, metadata, and values are
tested against the approved 2026.7 dependency.

`jwst_gtvt` supplies current visibility calculations, but legacy
`visibilityPA.py` still opens the packaged short ephemeris. That file cannot be
removed without first retiring or migrating the legacy public path.

## 10. Data retained unchanged

Retain the Fortney and generic databases unchanged for v2026.7, with new
manifests recording their verified payload checksums and unresolved builder
provenance. Retain ACES and ATLAS9 raw intensity grids as optional,
independently versioned limb-darkening components. Their scientific semantics
do not depend on Pandeia 2026.7 merely because the release number does.

Retain small active packaged auxiliary inputs and test fixtures until focused
PRs either document, regenerate, or prove them obsolete.

## 11. Data requiring regeneration or revalidation

* Regenerate the packaged groups table against the approved Pandeia 2026.7
  environment after recovering/writing a reproducible builder.
* Regenerate packaged JWST throughputs or prove runtime dependency equivalence.
* Regenerate complete, path-neutral ModelGrid indexes/caches if they will be
  distributed.
* Regenerate supported contamination traces only after all blockers clear.
* Regenerate order-0 only after scientific approval; do not infer approval from
  the existence of current operational files.

## 12. Newly required products

Current main newly exposes NIRCam DHS F322W2 and F444W contamination. Runtime
maps these apertures to parent `NRCA5_GRISM256_F322W2` and
`NRCA5_GRISM256_F444W` temperature-grid FITS and places the template into ten
DHS stripes. The local candidates occupy 26,812,800 and 36,792,000 bytes,
respectively, but the Box contamination ZIP lacks both directories. Their
release output size is provisional until regenerated.

The target-specific HDF5 cache is also new since the previous release, but it
is a runtime product, not a distributable data archive. The packaged
`contam_precompute_targets.txt` is a new generation input.

## 13. Obsolete products

Concrete removal candidates are limited to families with no current consumer:

* future external `groups_integrations.json`, because all current website and
  test callers pass the packaged table; only stale tutorial text reads it;
* the `all.zip` empty-directory placeholder;
* `__MACOSX`, `.DS_Store`, and AppleDouble files;
* packaged 2MASS J/H/K curves and the long 2018 ephemeris, for which exact-name
  and data-path searches found no current route, public API, worker, docs
  example, or test reader.

The many legacy contamination modes are **not** approved for removal. The web
marks several modes visibility-only, but the broad public `field_simulation`
aperture path remains. A support-contract PR must resolve that ambiguity first.

## 14. Mutable runtime state

`exoctk_log.db` is created from an in-code schema and populated by web form
submissions. Published log DBs differ. Future deployments should create or
migrate a runtime database outside immutable data components.

Contamination cache files are opened in append mode, receive new target groups
on misses, and depend on live target/catalog queries, date/visibility,
templates, pySIAF, code, and source classification. Current groups record only
name, coordinates, a filled flag, good PAs, target traces, sparse contaminant
planes, and plane indices. They lack schema, code, dependency, query, catalog,
epoch, and template versions.

## 15. Contamination-product dependencies and blockers

The web form exposes eight modes. Full contamination is enabled only for
SOSS SUBSTRIP256, SOSS SUBSTRIP96, and NIRCam DHS F322W2/F444W. NIRCam grism
F322W2/F444W, MIRI LRS, and NIRSpec are labeled visibility-only. Tests assert
this split. Generator, consumer, and old-data mode sets do not coincide.

Current trace scaling includes:

* SOSS trace values below 1 zeroed, then all orders multiplied by 1.5;
* contaminant traces multiplied by Gaia-derived relative flux and per-order
  empirical scales;
* NIRISS order-0 templates multiplied by `1.5e3`, then relative flux and the
  order-0 empirical scale;
* galaxy traces represented by masks rather than a dependency stellar SED.

Release blockers are mandatory and cumulative:

1. PR #725 must merge/finalize Gaia contaminant classification.
2. PR #726 must merge/finalize endpoint-aware Gaia TAP behavior.
3. Order-0 absolute calibration must receive explicit scientific sign-off.
4. Trace/cache manifests must identify code commit, Pandeia engine and refdata,
   pySIAF, SED/normalization, schema, template hashes, catalog/query behavior,
   and generation epoch.
5. Cache invalidation/migration must be defined before deployment regeneration.

No G-band/flux cutoff is proposed as a substitute for correct modeling.

## 16. Unresolved provenance questions

* Where are the original Fortney database builder and exact source models?
* Who built the generic grid, with what code, opacity inputs, and licenses?
* Which Pandeia engine/refdata built each groups and contamination product?
* Where were the local NIRCam F322W2/F444W templates published?
* Which original Kurucz files produced the ATLAS9 conversion?
* What source/version produced `predicted_gaia_colour.txt` and the SOSS
  wavelength mapping?
* Is local `output.zip` byte-identical to the current Box `modelgrid.zip`?
* Which cache set and external archive set are mounted in production?

## 17. Unresolved scientific questions

* What observational calibration justifies the order-0 template morphology
  and `1.5e3` amplitude?
* Are current empirical per-order/SOSS 1.5 scales valid for regenerated
  Pandeia 2026.7 products?
* Are local DHS parent traces scientifically valid for both filters and ten
  stripe placements?
* What numerical tolerances establish throughput equivalence?
* Can ACES/ATLAS9 be compressed or served on demand without changing
  interpolation and limb-darkening results?
* Which non-web contamination apertures remain supported public science APIs?

## 18. Proposed component layout

This is a manifest design recommendation, not an implemented installer:

* `tests`: small checked-in model-grid and contamination fixtures only.
* `core`: packaged small runtime/reference tables and website assets; groups
  may be a named subcomponent even if physically packaged.
* `website-smoke`: composition of core, Fortney, generic, ATLAS9, and the one
  SOSS trace set actually exercised by website CI. It should not download
  unrelated modes or ACES merely to run smoke tests.
* `forward-models-fortney` and `forward-models-generic`: separate because their
  provenance and sizes differ.
* `limb-darkening-atlas9` and `limb-darkening-aces`: separate; ACES is the
  dominant installation cost.
* `contamination-soss` and `contamination-nircam-dhs`: separate generated
  products with a shared schema/provenance contract.
* `website-current`: a composition containing every data-backed choice the
  deployed site actually exposes. If ACES remains selectable, it must include
  ACES even though smoke CI need not.
* `all`: a composition only, never an independent duplicate archive.

Logging and target cache directories are runtime volumes, not components.
Dependency-owned SED/refdata are referenced by dependency versions, not copied.

## 19. Estimated installation sizes

Known clean old-data sizes are 23.28 GB total, dominated by ACES (16.82 GB),
contamination (4.03 GB mixed), and generic (2.31 GB). The proposed component
sizes before regeneration are approximately:

| Component | Known size |
|---|---:|
| Tests/package-data baseline | 24.0 MB under `exoctk/data`; 27.6 MB including the contamination HDF5 fixture |
| Groups packaged table | 0.253 MB |
| Fortney | 46.50 MB |
| Generic | 2.313 GB |
| ATLAS9 local clean | 64.68 MB |
| ACES local clean | 16.824 GB |
| Current SOSS SUBSTRIP256 candidates | 1.260 GB |
| Current DHS parent candidates | 63.60 MB |
| Order 0 | 0.061 MB |
| Website-smoke composition | Roughly 3.7 GB using current candidate sizes |
| Website-current/all | Unknown until regeneration and support decisions; at least ~20.5 GB if both limb grids and current route families remain local |

Future contamination and groups sizes must not be presented as final until
actual validated products exist.

## 20. Recommended focused follow-up PR sequence

1. Merge this audit/inventory/tooling branch.
2. Add manifest schema, checksums, safe installer, component composition, and
   archive path-safety tests without changing scientific products.
3. Fix wheel/sdist package-data inclusion; separate log/cache runtime
   initialization; remove only proven unused packaged files and stale groups
   documentation.
4. Recover/write the groups builder, record Pandeia provenance, decide direct
   calls versus packaged table, regenerate if retained, and validate all modes.
5. Decide ACES/ATLAS9 distribution, cache/index schema, compression/on-demand
   feasibility, and provenance manifests.
6. Preserve Fortney/generic unchanged behind separate manifests while
   investigating their builders and licenses/provenance.
7. Add contamination cache schema, provenance, invalidation, and migration;
   do not regenerate target caches yet.
8. Merge/finalize PRs #725 and #726 independently.
9. Complete and approve order-0 scientific calibration in a focused PR/review.
10. Regenerate SOSS and DHS contamination products with manifests; validate
    mode lists, scaling, shapes, source rendering, and on-sky behavior.
11. Replace anonymous Box CI commands only after real component archives,
    checksums, and publication locations exist. Remove the duplicate and empty
    placeholder then.
12. Publish final v2026.7 components and release documentation.

## 21. Exact maintainer commands for unresolved checks

### Directly verify the large Box modelgrid object

This downloads 10.52 GB, requires about 11 GB temporary space for inspection
without extraction (about 29 GB if also extracted), and may take 10–90 minutes:

```bash
audit_tmp=$(mktemp -d /tmp/exoctk-modelgrid-audit.XXXXXX)
python scripts/audit_exoctk_data_archives.py \
  --url https://stsci.box.com/shared/static/vcn1ff9nzg2qx3pnkmsz17dujjdzlsh1 \
  --download --allow-large --output-dir "$audit_tmp" \
  --json "$audit_tmp/modelgrid-report.json"
shasum -a 256 "$audit_tmp/modelgrid.zip" /Users/tbell/exoctk_data/output.zip
```

Do not extract into `/Users/tbell/exoctk_data`.

### Inspect the remaining multi-GB unversioned STScI archives

The contamination, generic, and two ACES archives total about 13.69 GB
compressed. Metadata-only inspection is normally sufficient because Box/local
counterparts were inspected. If exact tar membership is required, budget at
least 14 GB temporary space without extraction and 30+ GB with extraction:

```bash
audit_tmp=$(mktemp -d /tmp/exoctk-stsci-audit.XXXXXX)
python scripts/audit_exoctk_data_archives.py \
  --url https://data.science.stsci.edu/redirect/JWST/ExoCTK/compressed/exoctk_contam.tar.gz \
  --url https://data.science.stsci.edu/redirect/JWST/ExoCTK/compressed/generic.tar.gz \
  --url https://data.science.stsci.edu/redirect/JWST/ExoCTK/compressed/modelgrid_ACES_1.tar.gz \
  --url https://data.science.stsci.edu/redirect/JWST/ExoCTK/compressed/modelgrid_ACES_2.tar.gz \
  --download --allow-large --output-dir "$audit_tmp" \
  --json "$audit_tmp/stsci-large-report.json"
```

### Record the exact deployment baseline

Run on the production host with read-only container commands and substitute the
actual Compose project/container names:

```bash
docker compose ps
docker compose images
docker compose exec exoctk python -c \
  'import exoctk, importlib.metadata as m; print(exoctk.__version__); print(m.version("pandeia.engine"))'
docker compose exec exoctk sh -c \
  'printf "%s\n" "$EXOCTK_DATA" "$EXOCTK_CONTAM_CACHE" "$pandeia_refdata"'
docker compose exec exoctk git -C /exoctk/exoctk rev-parse HEAD
find "$EXOCTK_DATA" -type f -exec shasum -a 256 {} + | sort
find "$EXOCTK_CONTAM_CACHE" -type f -exec shasum -a 256 {} + | sort
```

Do not infer the SHA if the deployed image has no `.git`; record that it is
unavailable and add immutable build labels in the manifest PR.

### Reproduce package-distribution inclusion check

```bash
audit_tmp=$(mktemp -d /tmp/exoctk-package-audit.XXXXXX)
git archive HEAD | tar -x -C "$audit_tmp"
python -m pip wheel --no-deps --no-build-isolation "$audit_tmp" -w "$audit_tmp/dist"
unzip -l "$audit_tmp"/dist/*.whl | rg 'exoctk/(data|exoctk_app/templates|exoctk_app/static)/'
(cd "$audit_tmp" && python setup.py sdist --dist-dir "$audit_tmp/dist")
tar -tzf "$audit_tmp"/dist/exoctk-*.tar.gz | rg '/exoctk/(data|exoctk_app/templates|exoctk_app/static)/'
```

### Validate a future groups table

The builder must first be recovered or written. Once it exists, the maintainer
should run it in a pinned Pandeia 2026.7 environment into a new temporary
directory, never over the tracked file, then run schema/value regression:

```bash
audit_tmp=$(mktemp -d /tmp/exoctk-groups-audit.XXXXXX)
python path/to/recovered_groups_builder.py --output "$audit_tmp/groups_integrations_input_data.json"
python -m json.tool "$audit_tmp/groups_integrations_input_data.json" >/dev/null
python -m pytest exoctk/tests/test_groups_integrations.py -q
```

The follow-up PR must add explicit engine/refdata version and generator tests;
no builder path is invented by this audit.

### Validate future throughput dependency equivalence

```bash
python -m pytest exoctk/tests/test_throughputs.py exoctk/tests/test_limb_darkening.py -q
```

Before deleting tables, add a new parameterized test that maps every packaged
JWST filename to a Pandeia instrument configuration and compares dimensions,
micron wavelength coverage, dimensionless PCE, metadata, and reviewed numeric
tolerance.

## Audit validation expectations

The audit branch validation includes JSON parsing, synthetic ZIP/tar helper
tests, metadata-only URL discovery/probing, Flake8 for changed Python, full
diff whitespace checks, package artifact inspection, old-tree before/after
metadata hashing, and review against `origin/main`. The final branch must
contain only this report, the JSON inventory, the read-only script, and its
small tests; no downloaded archive or generated scientific product belongs in
Git.
