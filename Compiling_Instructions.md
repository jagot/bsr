# Compilation Instructions for Various Systems
# Note: These assume everything is up to date, cmake may need updated as there is a minimum version requirement

## Anvil

```bash
ssh -XY x-ampuser@anvil.rcac.purdue.edu

# set up directories/modules/compilers
cd apps
mkdir bsr && cd bsr
git clone https://github.com/jagot/bsr
cd bsr
git checkout fix-compile
mkdir build && cd build
module load intel/19.1.3.304
export FC=$(which mpiifort)
export CC=$(which mpiicc)
export CXX=$(which mpiicpc)

# compile
cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_Fortran_COMPILER=mpiifort ../src/
make

# information about getting .sh files in /build/bin
# copy .sh files to local spot on machine
scp ampuser@login.expanse.sdsc.edu:/home/ampuser/apps/BSR/bsr/build/bin/*.sh /local/spot/on/machine

# copy .sh files to the super computer used in the build/bin directory
scp /local/spot/on/machine/*.sh /the/super/computer/used/:.../build/bin
```

## Bridges2
```bash
ssh ampuser@bridges2.psc.edu

# set up directories/modules/compilers
cd apps
mkdir bsr && cd bsr
git clone https://github.com/jagot/bsr
cd bsr
git checkout fix-compile
mkdir build && cd build
module load intel/2021.3.0
module load intel-mpi
export FC=$(which ifort)
export CC=$(which icc)
export CXX=$(which icpc)

# compile
cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_Fortran_COMPILER=mpiifort ../src/
make

# information about getting .sh files in /build/bin
# copy .sh files to local spot on machine
scp ampuser@login.expanse.sdsc.edu:/home/ampuser/apps/BSR/bsr/build/bin/*.sh /local/spot/on/machine

# copy .sh files to the super computer used in the build/bin directory
scp /local/spot/on/machine/*.sh /the/super/computer/used/:.../build/bin
```

## Expanse
```bash
ssh -X ampuser@login.expanse.sdsc.edu

# set up directories/modules/compilers
cd apps
mkdir BSR_updated && cd BSR_updated
git clone https://github.com/jagot/bsr
cd bsr
git checkout fix-compile
mkdir build && cd build
module load cpu/0.15.4
module load intel/19.1.1.217
module load intel-mkl/2019.1.144

# compile
cmake ../src/
make

# information about getting .sh files in /build/bin
cp ~/apps/BSR/bsr/build/bin/*.sh ~/apps/BSR_updatead/bsr/build/bin
```

## Frontera
```bash
ssh ampuser@frontera.tacc.utexas.edu

# set up directories/modules/compilers
cd apps
mkdir bsr && cd bsr
git clone https://github.com/jagot/bsr
cd bsr
git checkout fix-compile
mkdir build && cd build
module load intel
export FC=$(which mpif90)
export CC=$(which mpicc)
export CXX=$(which mpicxx)

# compile
cmake ../src/
make

# information about getting .sh files in /build/bin
# copy .sh files to local spot on machine
scp ampuser@login.expanse.sdsc.edu:/home/ampuser/apps/BSR/bsr/build/bin/*.sh /local/spot/on/machine

# copy .sh files to the super computer used in the build/bin directory
scp /local/spot/on/machine/*.sh /the/super/computer/used/:.../build/bin
```

## JetStream2
**In Progress**

## Ookami
```bash
ssh amos@login.ookami.stonybrook.edu

# set up directories/modules/compilers
cd apps2
mkdir bsr && cd bsr
git clone https://github.com/jagot/bsr
cd bsr
git checkout fix-compile
mkdir build && cd build
ssh fj-epyc
module load intel/compiler
module load intel/mpi
export FC=$(which ifort)
export CC=$(which icc)
export CXX=$(which icpc)

# compile
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ../src/
make

# information about getting .sh files in /build/bin
# copy .sh files to local spot on machine
scp ampuser@login.expanse.sdsc.edu:/home/ampuser/apps/BSR/bsr/build/bin/*.sh /local/spot/on/machine

# copy .sh files to the super computer used in the build/bin directory
scp /local/spot/on/machine/*.sh /the/super/computer/used/:.../build/bin
```
