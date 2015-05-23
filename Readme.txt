(1) Magudi uses a build system based on cmake. You will need to install it first. apt-get/brew/port install cmake should work on Linux/Mac.
(2) Create a separate directory for building that is not inside the source root. This is called an out-of-source build and is the preferred method to build my code.
(3) Run cmake <path to source directory>
(4) Add -DCMAKE_BUILD_TYPE=release to build with optimization enabled. By default, it builds in DEBUG mode (the DEBUG mode is significantly slower, more so than plascomcm since I use plenty of assertions to help with debugging).
(5) Now run make and this should build the executable.
(6) "make clean" cleans the built objects but keeps the cmake-generated build scripts including the Makefile. "make purge" will leave the directory in the state it was before you did anything with cmake removing all make/cmake-related files.
