
mkdir build-conda
cd build-conda

cmake .. ^
-D CMAKE_BUILD_TYPE:STRING="Release" ^
-D BUILD_TESTS:BOOL=FALSE ^
-D BUILD_PYTHON_BINDING=TRUE ^
-D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
-G "Ninja" ^
-Wno-dev

ninja install

cd ..
