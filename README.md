# Stage-wise Extreme Atomic Proximity in HPacker Side-chain Reconstruction

**Version 1.1.0 — 28 October 2026**  
**Version-specific Zenodo DOI (reserved for this release): [10.5281/zenodo.23717073](https://doi.org/10.5281/zenodo.23717073)**  
**Historical version 1.0.0 DOI:** [10.5281/zenodo.23616075](https://doi.org/10.5281/zenodo.23616075)

Author: Sho Miyazaki — Independent researcher, Tokyo, Japan. This research record reports 200 BioEmu WT source frames across six archived HPacker side-chain-reconstruction stages, and a strictly defined nonadjacent heavy-atom extreme-proximity threshold of 1.6 Å.

The E1 observation is 147/200 strict same-atom recontacts and 43/200 PHE6–PHE55 residue-pair ordered recontacts. Selected E2 local interventions/NN saved-array diagnostics, E3 four restrained short OpenMM continuations, and E4 a restricted protein heavy-atom LJ proxy are auxiliary records, **not** generalized causal, MD-success or free-energy claims.

## Manuscript, data and reproducibility

- [Technical Note v1.1.0](manuscript/Technical_Note_v1.1.0.pdf)
- [Supplementary Information v1.1.0](supplement/Supplementary_Information_v1.1.0.pdf)
- [Full version 1.1.0 Zenodo record](https://doi.org/10.5281/zenodo.23717073): all 1,200 stage PDBs and detailed E2–E4 evidence (once published)
- [Reproducibility instructions](REPRODUCIBILITY.md)
- [Provenance](PROVENANCE.md) and [rights scope](LICENSES/LICENSE_SCOPE.md)
- [Changes from v1.0.0](CHANGELOG.md)

The GitHub package is a lightweight entry point. Its historic v1.0.0 manuscript, SI, code and figures remain available at their preexisting paths. The v1.1.0 manuscript and SI have versioned filenames.

The full Zenodo package supports read-only verification with Python 3.10+ and NumPy:

```sh
python code/verify_package_hashes.py .
python code/verify_e1_from_pdb.py --root . --out verification/recheck
python code/verify_supporting_records.py --root . --out verification/recheck
```

E2 evidence is saved-array arithmetic and a single previously archived forward receipt; it is not ten independent full forward replays. E3 consists of selected, restrained OpenCL/single-precision continuations. E4 is a restricted LJ proxy rather than a thermodynamic free energy.

## Rights and acknowledgment

File-scoped rights, BioEmu and H-Packer upstream attribution, and third-party notices are documented in [NOTICE.md](NOTICE.md), [LICENSES/LICENSE_SCOPE.md](LICENSES/LICENSE_SCOPE.md), and the rights manifest. Author-owned text/figures use CC BY 4.0; eligible original code uses MIT; BioEmu coordinate lineage and third-party code retain applicable upstream terms. AI assistance and author responsibility are documented in the Technical Note.

**Citation:** cite the v1.1.0 version DOI 10.5281/zenodo.23717073 for these particular files, not the historical v1.0.0 version DOI.
