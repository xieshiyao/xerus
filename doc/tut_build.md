# Building Xerus                                 {#tutorial_build}

## Obtaining Xerus

You can get the source-code of the `xerus` library via [git](https://git.hemio.de/xerus/xerus/tree/master) or as an archived download under the same link. For example to 
clone the repository with the latest stable version under linux simply type
~~~.sh
git clone https://git.hemio.de/xerus/xerus.git
~~~

After downloading the source it is necessary to set a number of options that are somewhat individual to your system and needs. All following build steps assume that these
options are set in the `config.mk` file. Reasonable default values can be found in the `config.mk.default` file.

In particular the optimization level is interesting and the paths or standard flags for the required libraries need to be set. This way the `xerus` library can be compiled with any
`blas` compatible library for the matrix operations. For more details see the description in the `config.mk.default`
~~~.sh
cp config.mk.default config.mk
nano config.mk
~~~


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
Note in particular, that all 94 out of 94 tests were passed. Should this not be the case please file a bug report with as many details as you 
can in our [issuetracker](https://git.hemio.de/xerus/xerus/issues).


## Building the Library


## Compiling own Applications Using Xerus
