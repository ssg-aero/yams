*Build*

``` shell
rattler-build build --recipe yams/recipe/recipe.yaml -c ssg-aero -c conda-forge -m yams/recipe/variant_config.yaml
rattler-build upload anaconda output/linux-64/yams-*-py* --owner=ssg-aero -a=XXXXXXXXXXXXXXXXXXXX
```