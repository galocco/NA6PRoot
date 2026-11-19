#!/bin/bash
rm -rf build
rm -rf install
export ZLIB_ROOT=$HOME/alice/sw/ubuntu2404_x86-64/zlib/v1.3.1-4
export CMAKE_PREFIX_PATH="$ZLIB_ROOT:$CMAKE_PREFIX_PATH"
export XERCESC_ROOT="$HOME/alice/sw/ubuntu2404_x86-64/xercesc/Xerces-C_3_2_5-14"
export CMAKE_PREFIX_PATH="$XERCESC_ROOT:$CMAKE_PREFIX_PATH"
NA6PRootSrc=`pwd`;
BUILDDIR="$NA6PRootSrc/build"                                     # assuming that we want the build directory to be here ...
INSTALLDIR="$NA6PRootSrc/install"                                 # assuming that we want the installation directory to be here
mkdir -p $BUILDDIR; cd $BUILDDIR;
cmake -DCMAKE_INSTALL_PREFIX=$INSTALLDIR $NA6PRootSrc
make -j5 install    
