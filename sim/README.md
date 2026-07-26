<!-- doxy
\page refSIM sim
/doxy -->

# Simulation

Use executable
```
na6psim <options> --configKeyValues "<semicolon separated configurable params>"
```

### Layout definition
The detector layout is defined by the ConfigurableParam-based [NA6PLayoutParam](../base/include/NA6PLayoutParam.h) class. Any parameter defined in this class can be modified from the command line. For instance,
to set the `Z` of the 1-st VerTel station to `8.0` cm instead of its default value and then to move the whole VerTel in X and Z globally by `0.1` and `-0.2` cm respectively and change the 2nd target material to Iron (make sure that
corresponding material and medium are defined in the `createMaterials()` method) one can pass an option:
```
na6psim <options> --configKeyValues "layout.shiftVerTel[0]=0.1;layout.shiftVerTel[2]=-0.2;layout.posVerTelPlaneZ[0]=8.0;layout.medTarget[1]=Iron;<other settings>"
```
After the execution a corresponding `geometry.root` is generated in the working directory as well as `na6pLayout.ini` text file with the settings used for layout (the directory where this file is generated can be modified by `keyval.output_dir=<directory>` configurable).
Instead of passing long configKeyValues options string on the command line one can generate this `na6pLayout.ini` once (even with default options) and then edit the values of the `[layout]` block.
To run simulations with this modified layout (as well as with other configurable parames) one can use
```
na6psim <options> --load-ini <user-file.ini> --configKeyValues "<semicolon separated configurable params>"
```
Note that even in this case, one can still modify some parameters via `--configKeyValues ...` settings: the priority is given to the value provided from the command line, then loaded from the ini file (if any), and then to the class default setting.

Writing of the ini file can be disable by the option `--disable-write-ini`.

## Event-level parallel simulation

For event-level parallelism use the wrapper executable `na6psim_parallel`. It accepts the same command line options as `na6psim`, plus the wrapper option `--workers <N>`
which defines how many independent `na6psim` processes are run in parallel. The default is `--workers 1`.

Example:
```
na6psim_parallel --workers 5 -n 20 -g genDimuonBgEvent.C --configKeyValues "keyval.output_dir=sim_out"
```

The wrapper splits the requested event count into contiguous chunks (defined by chunk-specific `--event-offset`, `--skip-events`).
Each worker runs a separate VMC/TGeant4 process with its own temporary `keyval.output_dir`, `--event-offset`, `--skip-events`, and `--doDigitization false`.
After all workers finish, ROOT output files are merged into the final `keyval.output_dir`; `geometry.root` and `na6pLayout.ini` are copied from the first worker.
Unless the original command contains `--doDigitization false`, the wrapper then runs a final `na6psim --digitize-only`
step in the merged output directory, so digit MC labels use the final event numbering.
If `--digitize-only` is provided explicitly to `na6psim_parallel`, no workers are started and the wrapper runs a single `na6psim --digitize-only` command in the requested output directory.

Parallel mode requires an explicit non-negative `-n` / `--nevents` value, since the wrapper must split the event range before starting workers.

If the original command provides a non-negative random seed with `-r` / `--rnd-seed`, this seed is passed unchanged to the first worker and incremented by one for each following worker.
For example, `-r 100` with `--workers 3` runs workers with seeds `100`, `101`, and `102`.
If no or a negative seed is provided, each `na6psim` process generates its own seed from a hash of the current nanosecond time and the process ID.
File-backed generators such as `NA6PGenHepMC` support `--skip-events`;
generators which do not override event skipping will run different sequences thanks to modified (or time/process randomized) random seeds.

Before launching workers, `na6psim_parallel` prints the full list of worker commands, their log files, the temporary worker directory,
and the final merge/digitization directory. Use `--keep-worker-output` to preserve temporary worker directories for inspection after the merge.

## Available generators

Generators should be configured in root compilable macros with the generator producer function name being the same as the macro name.
The function signature should be `NA6PGenerator* <name>()`. The macro file name must be passed as `na6psim -g <macro_file_name>+ ...`

### GenBox

### GenParam

Allows to generate single type of particle with its transverse (pt or Mt) and longitudinal (Y or eta) distributions provided via string which are converted to TFormula.
The number of particles to generate per event is either provided explicitly (`setNTracks`) or estimated from provided `setdNdY` (or `setdNdEta`). This number is treated as
as Poisson mean in `setPoisson` was called.
See an example [test/genParamPi.C](../test/genParamPi.C).
The [test/genBgEvent.C](../test/genBgEvent.C) shows an example of CockTail generator production PbPb hadronic (pi, K, p) event for specific SPS energies.

### GenHepMC

Allows to import events in HEPMC format from external root file with `hepmc3_tree`. The generation can be done via wrapper macro [test/genHEPMC.C](../test/genHEPMC.C), e.g.
```
na6psim -n <N> -g $NA6PROOT_ROOT/share/test/genHEPMC.C+\(\"<input_file_name>\",false\)
```
In this example, events from the file `<input_file_name>` will be imported, discarding the intermediate particles already decayed by the generator.
The default behaviour is to add them as primaries (without further propagation). If the requested number of events `<N>` is negative or exceeds the number of events in the input file,
then it will be internally overridden to the number of available events. Otherwhise, only the first <N> events will be imported.

## User hooks

User can provide a macro with the function which is executed at each entry and exit of certain methods, e.g. `NA60PMC::AddParticles` (which allows to modify know particle table and their decays)
or `NA60PMC::selectTracksToSave` (which defines which particles will be saved). In the latter case the function e.g. may get the stack and enforce storing some particles.
The signature of the function must be `int <name>(int entryPoint, bool entering_or_exiting)`. The negative returned value signals a fatal error. E.g. macro `testHooks.C` below can be passed as
`na6psim -u testHooks.C+`:

```
#include <fairlogger/Logger.h>
#include <TVirtualMC.h>
#include "NA6PMCStack.h"
#include "NA6PMC.h"
#include "NA6PSimMisc.h"

int testHooks(int arg, bool inout)
{
  if (arg == UserHook::ADDParticles) {
    LOGP(info, "Calling user hook from {} of the AddParticles method", inout ? "entry":"exit");
  } else if (arg == UserHook::SelectParticles) {
    if (!inout) {
      return 0; // do nothing
    }
    LOGP(info, "Calling user hook from selectTracksToSave method entrance");
    auto mc = (NA6PMC*)TVirtualMC::GetMC();
    auto stack = mc->getMCStack();
    int ntr = stack->GetNtrack(), nPrim = stack->GetNprimary();
    for (int i = nPrim; i < ntr; i++) {
      auto* part = stack->GetParticle(i);
      int idMoth = part->GetFirstMother();
      int mothPdg = -1;
      while (idMoth >= 0) {
        auto* currMoth = stack->GetParticle(idMoth);
        int absPdg = std::abs(currMoth->GetPdgCode());
        if (absPdg == 11 || absPdg == 22 || absPdg == 211 || absPdg == 130 || absPdg == 321 || absPdg == 2212 || absPdg == 2112) {  // stop if a "stable" particle is found in the ancestors
          break;
        }
        if ((absPdg > 400 && absPdg < 600) || (absPdg > 4000 && absPdg < 6000) || absPdg == 100443 || absPdg == 20443) {
          part->SetBit(UserHook::KeepParticleBit); // force saving particle i (and its ancestors)
          break;
        }
        idMoth = currMoth->GetFirstMother();
      }
    }
  } else {
    LOGP(error, "Unknown hook ID {}", arg);
    return -1;
  }
  return 0;
}
```

## Vertex generation

By default the interaction vertex will be randomly generated on one of the subtargets according to their declared cross-sections and thicknesses, see the `Target` block of the [base/include/NA6PLayoutParam.h](../base/include/NA6PLayoutParam.h).
One can override this behaviour by providing a user macro `<uservtx>.C` (via `-V` of `--user-vertex` options) with the function having signature `void <uservtx(float& x, float& y, float& z)` which assignes to `x,y,z` the wanted vertex position. An example of such a macro is
[test/vtxInDump.C](../test/vtxInDump.C), which generates an interaction vertex in the main absorber plug. It can be called e.g. as

```
na6psim -n <N> -g <your_generator> -V $NA6PROOT_ROOT/share/test/vtxInDump.C+
```
