# conda config --set anaconda_upload yes

# to add the gag
git tag -a 2.rc1 -m 'first 2 release candidate'
git push origin 2.rc1

To build the package:
conda build . -c bioconda -c conda-forge

To upload to anaconda:
anaconda upload --user iossifovLab snakeobject...tar.bz2
NOTE: for some reason this fail on wigtop1!! :(
(base) yamrom@zermatt conda.recipe % 
