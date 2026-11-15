# Triplet bound-pair velocity-relaxation benchmark

Scientific code accompanying **Velocity Relaxation of Triplet Bound Pairs: Accuracy of a Local Coupling Approximation** (manuscript PHYSSCR-158810).

This repository compares a full electronic bound-band vertex with the complete second-order local coupling derived from the same one-dimensional two-electron Hamiltonian. The couplings share the reference dispersion and zero-transfer vertex without fitting. Both are evaluated with the same prescribed equilibrium coordinate bath and initial population at each interaction strength.

Repository: https://github.com/tomtaylor1992911-lang/PHYSSCR158810-triplet-pair-benchmark

## Quick reproduction

Use Python 3.10 or newer in an isolated environment:

```bash
python -m pip install -r requirements.txt
python reproduce.py
python verify_reproduction.py outputs/pilot
python plot_results.py outputs/pilot
```

The default performs only the two V/t0=6, N_K=96 pilot generators. It writes fresh arrays, numerical gates, response curves, comparison metrics and the signed momentum-transfer identity into `outputs/pilot`. An existing output directory is never overwritten; choose a new directory with `--output` to repeat a run. Default input/output paths are resolved from the script directory; an explicit relative `--output` is resolved from the current working directory.

```bash
python reproduce.py --suite main --output outputs/main_run
python reproduce.py --suite all --output outputs/all_run
python verify_reproduction.py outputs/all_run
python finite_ring.py --all --output outputs/finite_ring_all
```

`main` covers V/t0=5, 6, 8, 12 with N_K=96, 192, 384. `all` additionally includes the V/t0=6 amplitude-halving control at N_K=384 and alternative-cutoff control at N_K=192, 384: 30 nonzero generators in total. The numerical settings and all physical conditions are explicit in `parameters.json`. The finite-ring command uses genuine occupation-basis translation-block eigenstates, ring sizes L=24,36,48 and four diagnostic transitions; it does not substitute geometric states for eigensolves.

## Core implementation

| File | Purpose |
|---|---|
| `model.py` | Band, projected/local vertices, series check and momentum grid |
| `bath.py` | Positive thermal coordinate spectrum, direct quadrature and FFT memory diagnostic |
| `solver.py` | Column generator, equilibrium, graph/zero-mode checks, spectral response and full integral |
| `analysis.py` | Refined response difference and corrected exact-index signed-Q budget |
| `finite_ring.py` | Independent finite-ring electronic validation |
| `reproduce.py` | Portable runner for declared parameter suites |
| `verify_reproduction.py` | Compare fresh results with retained numeric references |
| `plot_results.py` | Draw response curves from fresh saved arrays on a common time axis |

The reported error comparison at N_K=384 gives time underestimates of approximately 15.28%, 10.52%, 5.86% and 2.59% at the four main interactions. `reference_comparisons.csv` contains all 15 reported paired comparisons, including both controls, and `reference_convergence.csv` contains the saved adjacent-grid comparisons. `reference_finite_ring.csv` preserves the 48 independent matrix-element comparisons. These are numeric reference outputs, not substitutes for running the supplied solver. The optional plotter draws diagnostic curves from fresh arrays; it does not claim to recreate the manuscript figure layout.

## Release validation

The portable entry points were executed with Python 3.13.14, NumPy 2.5.1, SciPy 1.18.1 and Matplotlib 3.11.1. A separate pilot and the full declared suite (30 generators, 15 paired comparisons) passed the saved-reference check. All eight checked metrics in every comparison matched the retained references exactly in this environment. The independent finite-ring command reproduced 72 eigensolves, 48 matrix elements and the two documented coarse-ring failures; all finest-ring points passed. `VALIDATION.json` records the test scope and tolerances. Floating-point last digits may differ with platform and BLAS.

## Scope and retained limitations

- Triplet pairing and nearest-neighbour attraction are inputs. No microscopic triplet-over-singlet selection mechanism is established.
- The coupling is scalar in a fixed spin sector. This is not finite-density superconductivity, spin/nematic relaxation or hydrodynamics.
- The bath is an imposed equilibrium spectrum. Its energies and temperature scale with the pair bandwidth across the four main cases. The four-point ordering is not a universal error law or a material-specific prediction.
- The population treatment is single-band, weak-coupling and second order. The cutoff excludes real single-quantum breakup at the retained order only. Numerical gates and memory diagnostics do not provide an exact non-Markovian error bound.
- The integrated time measures relaxation of a prepared velocity bias; pair number is conserved. It is not a pair-survival lifetime or superconducting stiffness.
- The response extremum search refines sign changes on a finite grid; it is not a rigorous enclosure of every possible extremum. The time integral itself uses the full decaying spectrum.
- The signed-Q budget is an algebraic identity of two complete generators, not an independent scattering-channel deletion experiment. Its individual bins have stronger grid dependence than the total time.
- Two coarse-ring matrix elements at V/t0=5, L=24 exceed the fixed absolute tolerance of 1e-9 t0. Their errors are about 1.851696e-9 and 1.582684e-9 t0. They remain recorded; larger rings pass without changing the tolerance.
- The full earlier bath-validation JSON was not completely snapshotted. Its historical hash and later reconciliation cannot reconstruct the missing original full text. This portable release does not claim to repair that provenance gap.

## Data and provenance

The complete scientific record is publicly available as the **Supplemental_Data_PHYSSCR158810_v2.zip** asset in [release data-v2](https://github.com/tomtaylor1992911-lang/PHYSSCR158810-triplet-pair-benchmark/releases/tag/data-v2). Download that named asset from the release page. Its size is 116,236,642 bytes and its SHA-256 is `0f0a8873648df09cda009924e4626d48b8409d7008b9e02be580a31661fff47f`. It contains the original numerical inputs, saved generators and eigensystems, response arrays, validation records, retained failed checks, and figure/table source data. The automatically generated GitHub source-code archives do not contain this data asset.

The data ZIP is preserved byte-for-byte, including 557 ZIP members and its internal manifest of 556 payload files. Its historical README describes the local-only status at packaging on 13 November 2026; the public release on 15 November 2026 supersedes that availability sentence. Historical paths and workflow records remain unchanged. No DOI has been assigned. The release tag points to the validated portable code commit `c6ee350459a4e8972024970e8ebd5a539b2c336f`.

`PROVENANCE.md` describes extraction and portability changes. The electronic and bath modules are byte-identical to the validated modules. The population kernel and final-v2 analysis retain their scientific algorithms; only historical workspace routing and local metadata are replaced by this portable runner. The portable code files omit local author paths. The complete data archive preserves historical provenance paths; it contains no manuscript or peer-review correspondence files.

No reuse license has been selected for this release. Public visibility alone is not an additional grant of copyright permissions. Dependency licenses remain those of their respective projects.
