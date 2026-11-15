# Release provenance

This portable code distribution is derived from the validated numerical implementation and final corrected analysis. Historical archive paths are intentionally not reproduced in this public note. The internal source-to-release audit retains those detailed mappings.

| Release file | Source role | Change |
|---|---|---|
| `model.py` | `model.py`; SHA-256 `56455daf5eecc507dc88f2dc46d47679368aad944a9755d9b7ce7bec413926de` | Byte-identical scientific module; no local paths or credentials. |
| `bath.py` | `bath.py`; SHA-256 `1f4c69aef6eb5e2ce56e40af1bf16e3148145fc9f92551b8f71c6aaa6b7b3aa5` | Byte-identical scientific module; no local paths or credentials. |
| `solver.py` | `execute.py`; SHA-256 `8e46937746d2b1ae730cf2b5451030592f0d0860b28ee87ee86d44ee6f6a3626` | Exact solve function except removed path-dependent bath loader and code_snapshot copying. Numerical gates, rates, eigensystem, populations, integral and outputs preserved. |
| `analysis.py` | `analyze_results.py`; SHA-256 `e114b7e896a39d325a659a3d1b72a66a4ea6b15a5c4f7f92430de807398b7b6d` | Retained final-v2 response extremum and exact-index signed-Q binning algorithms; inject saved arrays instead of workspace path lookup and omit local provenance paths. |
| `parameters.json` | `最小运行清单.json`; SHA-256 `6d0afb811626a4ca07a09896759b8bc212eed6c900ea1ee27981c5c4e25d9a6b` | Numerical inputs retained; absolute directory and obsolete design-stage status removed. Grid lists restricted to actually executed instances. |
| `reference_comparisons.csv` | `comparison_final_v2.csv`; SHA-256 `15e68e57471a7a9be3afa3f8091bed17b6d5f35783f5cb7cba49031d1c5d6b31` | All 15 final-v2 comparison rows and numerical strings retained; two workspace source-path columns removed. |
| `reference_finite_ring.csv` | `finite_ring_matrix_elements.csv`; SHA-256 `18a77cc4fcc59e575bc2a269b427e3ff5c3eafc16fa4ac874cff16f818f1888c` | Byte-identical numeric reference CSV with no absolute/local-account path. |
| `reference_convergence.csv` | `convergence_final_v2.csv`; SHA-256 `d3528cd8189369453b52f6384f03e26dddc8862add5d930aeb2a913905351989` | Byte-identical numeric reference CSV with no absolute/local-account path. |
| `finite_ring.py` | `run_fock_validation.py`; SHA-256 `ed0b281d24729245e627c41b661e302f763e0085e31ee6126ba010ef2bc9c94b` | Core basis, Hamiltonian, phase, eigensolver and matrix-element functions copied verbatim. Compact independent CLI removes local historical gate-file metadata and retains coarse failures. |

The full bath-validation snapshot limitation and coarse-ring failures are retained in the README. The portable CLI does not claim to execute the original historical preflight release-file chain; instead it runs the same numerical gates on every fresh population solve.

## Complete data archive

The original complete data ZIP is now publicly available as the named asset in [release data-v2](https://github.com/tomtaylor1992911-lang/PHYSSCR158810-triplet-pair-benchmark/releases/tag/data-v2). The ZIP bytes and internal manifest are unchanged. Its historical local-only availability statement is superseded by the release; the documented missing earlier bath-validation JSON remains a provenance limitation. See the release notes for the file size and SHA-256.
