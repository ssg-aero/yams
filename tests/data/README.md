# Input grids of the C++ tests

Meridional grids read by `tests_curvaturesolver`. They are outputs of the tbs
tests (`tbs/tests/tests_merid_mesh.cpp`, cases `for_test005` to `for_test007`,
inputs in `tbs/tests/in`), generated with tbs f7d1953:

| File | tbs test |
|---|---|
| test_005.vts | merid_mesh.for_test005 |
| test_006.vts | merid_mesh.for_test006 |
| test_007.vts | merid_mesh.for_test007 |

The other grids the tests refer to (`test_001*`, `test_002` to `test_004`,
`test_008`, `test_009`, `alpx001`, `turbine_grid`) came from the former tbslib
and cannot be regenerated: their inputs are private data that could not be found.
The test cases that need them are skipped (`yams_tests::data_missing`).
