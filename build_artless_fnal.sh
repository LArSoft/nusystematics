source /grid/fermiapp/products/dune/setup_dune.sh
setup cmake v3_9_0
setup gcc v6_4_0

# setup ROOT
source /cvmfs/larsoft.opensciencegrid.org/products/root/v6_10_08b/Linux64bit+2.6-2.12-e15-nu-prof/bin/thisroot.sh

setup genie        v2_12_8c   -q e15:prof
setup genie_xsec   v2_12_8    -q DefaultPlusMECWithNC
setup genie_phyopt v2_12_8    -q dkcharmtau
setup dk2nu        v01_05_01b -q e15:prof
setup geant4 v4_10_3_p01b -q e15:prof
setup ifdhc

# we need something called TBB for ROOT
export TBB_DIR=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1
export TBB_LIB=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/lib
export TBB_INC=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/include

mkdir build
cd build
cmake ../ -DUSEART=0 -DLIBXML2_LIB=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/lib/ -DLIBXML2_INC=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/include/libxml2 -DPYTHIA6=/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28i/Linux64bit+2.6-2.12-gcc640-prof/lib
make systematicstools
make
make install
