# NA6PRoot framework

## Scope

NA6P framework for detector geometry definition, VMC-based simulation (and simulated data export for external packages) and reconstruction.

## Installation (to be completed)

Currently, it is assumed that the needed third-party packages (see `find_package` directives in the top level [CMakeLists.txt](CMakeLists.txt) for the list of these packages)
are installed e.g. via [O2](https://github.com/AliceO2Group/AliceO2) framework and one has `ALIBUILD_WORK_DIR` environmnet variable pointing to the ALICE `<...>/sw` directory.

After the installation of these packages some paths have to be initialized, e.g. by running (in case of installation via O2, bash shell is assumed)
```
myarch=$(aliBuild architecture)
alienv load O2/latest
export WORK_DIR=$ALIBUILD_WORK_DIR
source $WORK_DIR/$myarch/O2/latest/etc/profile.d/init.sh
```

Go to the directory you want to install the package and run
```
git clone https://github.com/shahor02/NA6PRoot.git                # clone the source
cd NA6PRoot;
export NA6PRootSrc=`pwd`;
export NA6PBUILDDIR="$NA6PRootSrc/build"                          # assuming that we want the build directory to be here ...
export NA6PINSTALLDIR="$NA6PRootSrc/install"                      # assuming that we want the installation directory to be here
mkdir -p $NA6PBUILDDIR; cd $NA6PBUILDDIR;
cmake -DCMAKE_INSTALL_PREFIX=$NA6PINSTALLDIR $NA6PRootSrc
make -j5 install                                                  # build and install
```

To run some `NA6PRoot` executable, you need to initialize its paths (as well as those of the packages `NA6PRoot` depends on, see above),
which particularly sets the NA6PROOT_ROOT.

```
source $NA6PINSTALLDIR/init.sh
mkdir tst
cd tst
na6psim -n 5 -g $NA6PROOT_ROOT/share/test/genbox.C+
```


## Monte Carlo simulation: [sim module](sim)
<!-- doxy
* \subpage refSIM
/doxy -->

## Reconstruction parameter checks

GitHub Actions validates every `.ini` file under `params/`, including
subdirectories, using Python's standard library. Run the same check locally:

```bash
python3 test/ci/checkRecoParams.py
```

No build, ROOT, O2PDPSuite, CVMFS, Docker, or extra Python packages are required.
The validator in `test/ci/checkRecoParams.py` checks INI syntax, duplicate keys,
the `[reco]` section, case-sensitive parameter names, scalar types, and array
indices against the declarations in `rec/include/NA6PRecoParam.h`. Partial
configurations are allowed; omitted parameters retain their application defaults.
It fails if no `.ini` files are found or any file is invalid, and reports errors
in the workflow log. This static check does not run simulation or reconstruction,
validate physical parameter choices, or inspect output trees.

The full integration test remains available locally after loading the project's
dependencies:

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel 5
cmake --install build
source install/init.sh
ctest --test-dir build --output-on-failure -R '^na6prec_params$'
```

That test generates one event with a fixed random seed and runs `na6prec
--load-recoparam` for every parameter file. A compiled checker verifies one event
in each cluster, vertex, VT/MS track, and matched-track output tree; empty track
collections are allowed. Failed runs keep their logs and outputs under
`build/Testing/params/`; successful runs are removed.

## Reconstruction: [rec module](rec)
<!-- doxy
* \subpage refREC
/doxy -->

