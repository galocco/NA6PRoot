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

Detailed detector documentation:

- [Vertex Telescope geometry and data flow](docs/VerTelGeometryAndDataFlow.md)
- [MWPC module and modular muon spectrometer](docs/MWPCMuonSpectrometerGeometry.md)

## Reconstruction: [rec module](rec)
<!-- doxy
* \subpage refREC
/doxy -->
