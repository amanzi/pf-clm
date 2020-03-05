
ORIGIN_DIR=$PWD
mkdir build ; cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$ORIGIN_DIR \
  -DBUILD_SHARED_LIBS=ON \
  -DUNSTRUCTURED_HOST=ON
make -j 1 VERBOSE=1
make install
cd $ORIGIN_DIR

