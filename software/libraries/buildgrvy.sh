#!/bin/sh

export GRVY_VERSION=0.32.0
export BOOST_VERSION=1.55.0
echo $BOOST_DIR

TOPDIR=`pwd`
rm -rf $TOPDIR/grvy-$GRVY_VERSION-src
tar xvzf grvy-$GRVY_VERSION.tar.gz
mv grvy-$GRVY_VERSION grvy-$GRVY_VERSION-src
cd grvy-$GRVY_VERSION-src || exit 1

#export GRVY_DIR=/workspace/pclark/Libraries/grvy-$GRVY_VERSION-$COMPILER-$COMPILER_VERSION-boost-$BOOST_VERSION

./configure --prefix=$GRVY_DIR --with-boost=$BOOST_DIR --without-hdf5 2>&1 | tee configure.log
make -j 4 2>&1 | tee make.log
make check 2>&1 | tee check.log
make install
mv config.log configure.log $GRVY_DIR
mv make.log check.log $GRVY_DIR
cd ..
#rm -rf $TOPDIR/grvy-$GRVY_VERSION-src
