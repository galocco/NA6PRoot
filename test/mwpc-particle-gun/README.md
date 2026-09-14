# MWPC particle-gun smoke test

This is a deliberately small transport test for the detailed `NA6PMWPCChamber`
geometry. It is the NA6PRoot-native analogue of the earlier isolated-chamber
muon tests in `na60-dice-mwpc-simulation`.

The test does **not** modify `NA6PMuonSpecModular` and does not define a station
layout. It creates one complete MWPC chamber at the origin, shoots one controlled
forward muon per event through the centre, and lets Geant4 transport the particle
through the ROOT/TGeo geometry via Geant4VMC.

## What it checks

- the installed branch can link and instantiate `NA6PMWPCChamber`;
- the complete 22-part chamber is passed to ROOT/TGeo and Geant4VMC;
- ROOT reports no overlaps in the isolated geometry;
- the chamber gas is recognized as a sensitive NA6PRoot volume;
- the primary muon enters the gas and produces Geant4 transport steps there;
- gas energy deposition is recorded as a simple transport-level diagnostic.

This is **not** an MWPC signal-response simulation. There is no avalanche gain,
drift, wire response, threshold, electronics, inefficiency or reconstruction.

## Prerequisite

Build and install the current `feature/mwpc-dice` branch first. The test links to
the installed `simLib`, `baseLib` and `utilsLib`, so this also checks that the new
chamber code is part of a real NA6PRoot build rather than only a standalone stub.

With the usual O2/NA6PRoot environment loaded:

```bash
cd NA6PRoot
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="$PWD/../install" ..
make -j5 install
source ../install/init.sh
```

If your installation lives elsewhere, source that installation's `init.sh` so
`NA6PROOT_ROOT` points to it.

## Run

From the repository root:

```bash
bash test/mwpc-particle-gun/run.sh
```

Defaults are 100 mu- events at 20 GeV/c with seed 20260914.

Arguments are:

```text
run.sh [events] [momentum_GeV] [PDG] [seed]
```

Examples:

```bash
# 100 central 20 GeV/c mu-
bash test/mwpc-particle-gun/run.sh 100 20 13

# 50 central 4 GeV/c mu+
bash test/mwpc-particle-gun/run.sh 50 4 -13
```

Each invocation creates a new time-stamped directory under `runs/`.

## Output

`summary.txt` contains the overall PASS/FAIL result. The run passes only when the
primary muon enters the sensitive gas and has gas transport steps in every event.

`events.csv` contains one row per event:

```text
event,gas_steps,gas_entries,gas_exits,primary_steps,primary_entries,energy_deposit_keV
```

`geometry.root` is the exact isolated TGeo geometry given to Geant4VMC.
`geometry_summary.txt` records the gas dimensions/centre and full chamber envelope.
`geometry_checks.txt` records ROOT's overlap result. `MCKine.root` is the usual
NA6PRoot Monte Carlo stack output. Empty hit files from other registered NA6PRoot
modules may also appear because the standard `NA6PMC::init()` lifecycle is reused;
they are not part of this test.

## Code flow

```text
mwpc_gun.C
    -> one known forward muon/event

MWPCTestModule
    -> calls the real NA6PMWPCChamber builder
    -> gas gets a normal NA6PRoot sensitive-volume ID

MWPCParticleGunApplication
    -> creates an air world + one chamber only
    -> exports geometry.root
    -> hands the TGeo geometry to Geant4VMC with geomRoot

Geant4
    -> transports the muon
    -> NA6PMC::Stepping() routes sensitive gas steps to MWPCTestModule::stepManager()
    -> events.csv / summary.txt
```

The local `MWPCTestModule` is test scaffolding only. Production integration should
still happen through the appropriate muon-spectrometer module (currently
`NA6PMuonSpecModular`) once the station/chamber layout is ready.
