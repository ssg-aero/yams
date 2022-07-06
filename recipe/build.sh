mkdir build-conda
cd build-conda

cmake .. \
-D CMAKE_BUILD_TYPE:STRING="Release" \
-D BUILD_TESTS=OFF \
-D BUILD_CLI=ON \
-D BUILD_PYTHON_BINDINGS=ON \
-D CMAKE_INSTALL_PREFIX=$PREFIX \
-G "Ninja" \
-Wno-dev

ninja install

cd ../python
$PYTHON -m pip install . --no-deps -vv

cd ..