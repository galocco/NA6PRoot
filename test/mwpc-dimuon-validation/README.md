# MWPC dimuon validation

This test uses the existing `test/genDimuonBgEvent.C` parent generator for three signal channels:

- `Jpsi` (PDG 443)
- `Omega` (PDG 223)
- `Phi` (PDG 333)

The validation hook `mwpcDimuonHooks.C` forces each signal parent to decay to `mu+ mu-` with 100% branching fraction for this test only. It also marks all signal daughter muons for storage in `MCKine.root`, including daughters which do not hit a sensitive detector volume. No production decay table is changed in the source code.

The default validation run uses 40 GeV/nucleon beam energy, one signal parent per event, `0 < pT < 5 GeV/c`, `0 < y < 6`, no hadronic background, and 20,000 events per channel. This gives 20,000 generated parents and 40,000 generated signal daughter muons per channel if the forced two-body decay works as intended.

Run from the repository after sourcing the normal NA6PRoot environment and selecting the branch install:

```bash
bash test/mwpc-dimuon-validation/run_validation.sh 20000 4
```

The second argument is the number of independent `na6psim` worker processes. The script first performs a one-event serial compile/smoke step and then runs the three channels sequentially with `na6psim_parallel`.

Outputs are written below:

```text
test_runs/mwpc_dimuon/Jpsi/
test_runs/mwpc_dimuon/Omega/
test_runs/mwpc_dimuon/Phi/
```

For each channel, `analyzeMWPCDimuon.C` produces:

- `momenta_<channel>.png`: parent |p| and pT, daughter-muon |p| and pT;
- `coverage_MS0_<channel>.png` ... `coverage_MS5_<channel>.png`: mean number of distinct A/B/C/D stagger layers hit by an accepted signal daughter muon, binned at the nominal station plane;
- `validation_<channel>.root`: the underlying histograms.

For a coverage bin, the analysis uses the nearest MWPC gas hit of that signal muon and linearly projects its incoming position/momentum over the small A/B/C/D z offset to the nominal station z. The bin value is the mean number of distinct checkerboard layers hit. Bins with fewer than 10 accepted signal muons are masked. The default bin size is 10 cm.

The coverage observable is therefore transport-level and conditional on the signal muon producing at least one MWPC hit in that station. It is intended as a direct validation of the implemented stagger/overlap response. If an exact reproduction of an external offline definition of "incoming muon" is required later, a dedicated truth scoring plane can be added.
