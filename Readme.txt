Merged follows the Pimpl Idiom (see http://www.c2.com/cgi/wiki?PimplIdiom)

Installation directions:
(1) Magudi uses a build system based on cmake. You will need to install it first. apt-get/brew/port install cmake should work on Linux/Mac.
(2) Create a separate directory for building that is not inside the source root. This is called an out-of-source build and is the preferred method to build Magudi.
(3) Run cmake <path to source directory>
    * Add -DCMAKE_BUILD_TYPE=release to build with optimization enabled. By default, it builds in DEBUG mode (the DEBUG mode is significantly slower).
    * Add -DPLOT3D_FORMAT=extended to allow for i/o on larger grids (>100 M grid points)
(4) Now run make and this should build the executable.
(5) "make clean" cleans the built objects but keeps the cmake-generated build scripts including the Makefile. "make purge" will leave the directory in the state it was before you did anything with cmake removing all make/cmake-related files.
