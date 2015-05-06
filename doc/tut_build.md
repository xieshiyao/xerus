# Building Xerus                                 {#tutorial_build}

## Obtaining Xerus

You can get the source-code of the `xerus` library via [git](https://git.hemio.de/xerus/xerus/tree/master) or as an archiv from [here](https://git.hemio.de/xerus/xerus/repository/archive.zip?ref=master).
For example to clone the repository with the latest stable version under linux simply type
~~~.sh
git clone https://git.hemio.de/xerus/xerus.git
~~~

## Make sure all dependencys are present
Xerus depends on several well established libraries, usually pre-installed or available through the standard package managers. In particular these `lapacke` and `cblas` and their repsective dependencies,
as well as `CXSparse` and `bfd` and their dependencies. Also Note that at at least version 4.8 of the `GCC` is required, as this is the first version to offer support for all `C++11` functuanality used in xerus. 
Make sure that all these are installed on your system before proceeding.



## Building the Library

After downloading the source it is necessary to set a number of options that are somewhat individual to your system and needs. All following build steps assume that these
options are set in the `config.mk` file. Reasonable default values can be found in the `config.mk.default` file.

In particular the optimization level is interesting and the paths or standard flags for the required libraries need to be set. This way the `xerus` library can be compiled with any
`blas` compatible library for the matrix operations. For more details see the description in the `config.mk.default`
~~~.sh
cp config.mk.default config.mk
nano config.mk
~~~

After setting the required options simply call make
~~~.sh
make all -j4
~~~
this creates the shared library object libxerus.so as well as the static library object libxerus.a. 

## Making Sure Everything Works as Intended
The sources include a number of unit tests that ensure that everything works as intended. To perform them simply input
~~~.sh
make test -j4
~~~
(where `-j4` allows make to use up to 4 threads for the compilation). If all options were set correctly in the `config.mk` it should compile a test executable and launch it.
The output of this executable should then list a number of passed tests and end with
~~~.sh
|
|                        total summary 94 of 94 passed                     
-------------------------------------------------------------------------------
|
|                       Total time elapsed: 5177.275 ms                    
-------------------------------------------------------------------------------
~~~
Note in particular, that all tests were passed. Should this not be the case please file a bug report with as many details as you 
can in our [issuetracker](https://git.hemio.de/xerus/xerus/issues).

## Installing xerus

If you want to install xerus on your system to the path given in `config.mk` simply call (as root is necessary)
~~~.sh
make install
~~~


## Compiling own Applications Using Xerus

If xerus is properly installed on your system, compiling your own applications using xerus is as simple as using any library. Just include `-lxerus` in your linker call and make sure to use
`-std=c++11` or `-std=c++14` in all compilation units that include `xerus.h` and everything should work.

If you want to use the static version of xerus you also have to include all libraries xerus depends opon in your linker call. In particular these are lapacke (-llapacke), 
lapack (-llapack), cblas (-lcblas) blas (-lblas), CXSparse (-lcxsparse), bfd (-lbfd). On some old system one has to manually include the dependencys of bfd (-liberty -lz -ldl).

You can test this by trying to compile the tutorial file
~~~.sh
g++ -std=c++11 -I misc tutorial.cpp libxerus.a -llapacke -llapack -lcblas -lblas -lcxsparse -lbfd -liberty -lz -ldl
~~~

