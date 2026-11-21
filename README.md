# shastry-sutherland-ed

MATLAB code for the exact diagonalisation (ED) of the spin-1/2 Heisenberg model on a finite cluster of 12 orthogonal dimers (24 spins) of the Shastry-Sutherland lattice, with the calculation of the corresponding magnetic neutron-scattering intensity.

The code accompanies the article:

> about to be submitted

## Model

```
H =  J  * sum_{dimers (i,j)}  S_i . S_j      (intra-dimer bonds)
   + J' * sum_{bonds  (i,j)}  S_i . S_j      (inter-dimer bonds)
```

The default parameters are `J = 7.023` and `J'/J = 0.63`, with the lattice constant `a = 8.985 Angstrom` (SrCu2(BO3)2 geometry). The cluster has open boundaries: 24 spins, 12 dimers, 32 inter-dimer bonds.

## Contents

| File | Description |
|---|---|
| `ED_SS_full_hilbert_space.m` | ED in the full Hilbert space (2^24 states). Includes a tiny Zeeman field so that every eigenstate has a definite S^z. Also computes the neutron cross section from the full ED wavefunctions. Needs a lot of memory. |
| `ED_SS_Sz_sector.m` | Memory-efficient ED in a fixed total-S^z sector. Includes the dimer-basis decomposition of a chosen eigenstate and the one-triplet-approximation cross section. |

## Requirements

- MATLAB (tested with: <MATLAB version>). No additional toolbox is needed.
- Memory:
  - `ED_SS_full_hilbert_space.m`: several tens of GB of RAM (the 72 sparse single-spin operators alone need about 30 GB).
  - `ED_SS_Sz_sector.m`: a few GB for N = 24 (sector dimension 2,496,144 for S^z = 1).

## Usage

Both scripts are organised in sections (`%%`). The parameters are at the top of each script.

1. Open the script in MATLAB.
2. Adjust the parameters in the first section:
   - `ED_SS_full_hilbert_space.m`: `J`, `J_prime`, `B`, `n_states` (number of eigenstates), `target_state` (state analysed below).
   - `ED_SS_Sz_sector.m`: `N` (8, 24 or 32), `Sz_target` (sector), `num_eig`, `target_state`, cutoffs of the dimer-basis decomposition, `save_tables`, `output_dir`.
3. Run the script. The spectrum is printed in the command window; figures are produced at the end of each analysis section.

The geometry (spin positions) is defined for the 24-spin cluster only. For other values of `N` in `ED_SS_Sz_sector.m` the script stops after the diagonalisation.

## What is computed

- Lowest eigenvalues and eigenstates (energy relative to the ground state, total S^z).
- Spin expectation values of a chosen eigenstate, plotted on the cluster.
- Decomposition of the eigenstate in the dimer basis (singlet `s`, triplets `t+`, `t0`, `t-` on each dimer), with the dominant components printed and optionally saved to Excel.
- Magnetic neutron-scattering intensity in the (h, k, 0) plane:
  - from the full ED wavefunctions of the ground state and the selected state (full-space script);
  - from the one-triplet components of the selected state, with the ground state taken as a product of dimer singlets (both scripts).

The intensity is `I(q) = sum_alpha |<f| sum_i exp(i q.r_i) S_i^alpha |0>|^2`. The polarisation factor, magnetic form factor and energy-dependent prefactors are not included.

## Conventions

- Spins are numbered 1..24; spins (2m-1, 2m) form dimer m.
- Dimer basis: 0 = singlet, 1 = t+, 2 = t0, 3 = t-.
- The two scripts index the basis states differently (kron ordering with digit 0 = up in the full-space script; bit encoding with bit 1 = up in the sector script) and use opposite signs for the singlet. Eigenvalues and intensities are not affected, only the sign of some dimer-basis amplitudes. Details are in the header of each script.
- `q = 2*pi/a * (h, k, 0)`, with `h`, `k` in reciprocal lattice units.


## Citation

If you use this code, please cite the article above (see also `CITATION.cff`).

## License

This code is released under the MIT license (see `LICENSE`).
