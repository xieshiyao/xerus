# Debugging

`xerus` comes with a number of functionalities to simplify debugging. The purpose of this page is to explain the most important ones.

## LOG

The library uses a macro of the form `LOG(level, msg)` to print messages to `cout`. It allows to use arbitrary log levels without any need to declare them beforehand
and to use typical piping syntax for the messages.
~~~.cpp
LOG(als_warning, "The ALS encountered a mishap " << variable_1 << " <= " << variable_2);
~~~
The following warning levels are predefined: `fatal`, `critical`, `error`, `warning`, `info`, `debug`. The default config file `config.mk.default` defines the preprocessor
variable `INFO_` which causes all but the `debug` level to be printed. Per default any other log-level is printed. This can be turned off by including 
`SET_LOGGING(als_warning, err::NO_LOGGING)` inside a commonly included header (like `include/xerus/tensorLogger.h`).

The `fatal` loglevel is special in that it will not just print the message but also throw an exception including the message itself and a callstack in the `.what()` string.

With the additional option `LOGGING += -D LOG_BUFFER_` in the `config.mk` file, any not printed log message will be stored in a circular buffer and will be dumped to a file
in the `errors/` subfolder (if it exists) on any `error`, `critical` or `fatal` log message. 


## REQUIRE

The `REQUIRE(condition, additional_msg)` macro replaces assertions in the `xerus` library. It uses above `LOG` functionality and can thus use piping style messages just as the `LOG`
macro. It is equivalent to the following definition
~~~.cpp
REQUIRE(condition, additional_msg) = if (!(condition)) { LOG(fatal, additional_msg); }
~~~
There is a large number of such checks in the library. All of them can be turned off by defining `DEBUG += -D DISABLE_RUNTIME_CHECKS_` in the `config.mk` file.


## UNIT_TEST

Compiling with `make test` creates an executable that includes all functions defined with the `UNIT_TEST(group_name, test_name, source_code...)` macro. This will only include the 
unit tests in the library but can easily be extended should you wish to use it for your own executable. The created executable can be launched manually to perform individual
tests (eg. `./XerusTest TTTensor:summation`) or reperform all of them (`./XerusTest all`).

With the additional `config.mk` option `DEBUG += -D TEST_COVERAGE_`, the executable will track which `REQUIRE` and `REQUIRE_TEST` macros were executed during the run to advice 
about function that need additional testing.


## Callstacks and XERUS_THROW

Unless `NO_FANCY_CALLSTACK = TRUE` was declared in the `config.mk` file, the function `get_call_stack()` will return a formatted string of the full stack trace including function name,
source file and line numbers. This stack will in particular be included in the exceptions thrown by `LOG(fatal, ...)` and `REQUIRE(...)` macros to simplify the detection of errors.

The information of this callstack is only available if the application was compiled with the `-g` flag and the linking requires the `binutils` packages `-lbfd -lz -ldl -liberty`.

The exceptions used by `xerus` have the additional capability to accept piped information that will be included in the `.what()` string. To include a callstack it is thus possible to 
simply write
~~~.cpp
XERUS_THROW(xerus::generic_error() << "callstack:\n" << xerus::get_call_stack());
~~~
The used macro will additionally include the source file and line as well as the function name in which the exception was thrown.

