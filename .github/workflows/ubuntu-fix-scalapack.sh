#!/bin/bash

# Try to fix the installation of ScaLAPACK on Ubuntu, cf.
# https://bugs.launchpad.net/ubuntu/+source/scalapack/+bug/1993843
# https://bugs.launchpad.net/ubuntu/+source/scalapack/+bug/1917534
# https://github.com/reference-scalapack/scalapack/issues/28

VERSION=$(dpkg-query --showformat='${Version}' --show libscalapack-openmpi-dev | cut -d '-' -f 1)
ARCH=$(arch)
src="/usr/lib/${ARCH}-linux-gnu/libscalapack-openmpi.so"
dst="/usr/lib/libscalapack-openmpi.so.${VERSION}"
echo "Trying to link ${src} -> ${dst}"
sudo ln -sv $src $dst
