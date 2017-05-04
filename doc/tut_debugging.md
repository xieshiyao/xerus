# Debugging

`xerus` comes with a number of functionalities to simplify debugging. The purpose of this page is to explain the most important ones.

## LOG

The library uses a macro of the form `XERUS_LOG(level, msg)` to print messages to `cout`. It allows to use arbitrary log levels without any need to declare them beforehand
and to use typical piping syntax for the messages.
~~~.cpp
XERUS_LOG(als_warning, "The ALS encountered a mishap " << variable_1 << " <= " << variable_2);
~~~
The following warning levels are predefined: `fatal`, `critical`, `error`, `warning`, `info`, `debug`. The default config file `config.mk.default` defines the preprocessor
variable `XERUS_INFO` which causes all but the `debug` level to be printed. Per default any other log-level is printed. This can be turned off by including 
`XERUS_SET_LOGGING(level, xerus::err::NO_LOGGING)` inside a commonly included header.

The `fatal` loglevel is special in that it will not just print the message but also throw an exception including the message itself and a callstack in the `.what()` string.

With the additional option `LOGGING += -D XERUS_LOG_BUFFER` in the `config.mk` file, any not printed log message will be stored in a circular buffer and will be dumped to a file
in the `errors/` subfolder (if it exists) on any `error`, `critical` or `fatal` log message. 


## REQUIRE

The `XERUS_REQUIRE(condition, additional_msg)` macro replaces assertions in the `xerus` library. It uses above `XERUS_LOG` functionality and can thus use piping style messages just as the `XERUS_LOG`
macro. It is equivalent to the following definition
~~~.cpp
XERUS_REQUIRE(condition, additional_msg) = if (!(condition)) { XERUS_LOG(fatal, additional_msg); }
~~~
There is a large number of such checks in the library. All of them can be turned off by defining `DEBUG += -D XERUS_DISABLE_RUNTIME_CHECKS` in the `config.mk` file.


## UNIT_TEST

Compiling with `make test` creates an executable that includes all functions defined within `xerus::UnitTest` objects.
~~~.cpp
static xerus::misc::UnitTest objectName("Group", "Name", [](){
    // source code of the test
    // likely includes (several) tests of the form TEST(condition);
    TEST(true);
});
~~~
This will only include the 
unit tests in the library but can easily be extended should you wish to use it for your own executable. The created executable can be launched manually to perform individual
tests (eg. `./XerusTest TTTensor:summation`) or reperform all of them (`./XerusTest all`).

With the additional `config.mk` option `DEBUG += -D XERUS_TEST_COVERAGE`, the executable will track which `XERUS_REQUIRE` and `XERUS_REQUIRE_TEST` macros were executed during the run to advice 
about function that need additional testing.


## Callstacks and XERUS_THROW

Unless `XERUS_NO_FANCY_CALLSTACK = TRUE` was declared in the `config.mk` file, the function `xerus::misc::get_call_stack()` will return a formatted string of the full stack trace including function name,
source file and line numbers. This stack will in particular be included in the exceptions thrown by `XERUS_LOG(fatal, ...)` and `XERUS_REQUIRE(...)` macros to simplify the detection of errors.

The information of this callstack is only available if the application was compiled with the `-g` flag and the linking requires the `binutils` packages `-lbfd -lz -ldl -liberty`.

The exceptions used by `xerus` have the additional capability to accept piped information that will be included in the `.what()` string. To include a callstack it is thus possible to 
simply write
~~~.cpp
XERUS_THROW(xerus::misc::generic_error() << "callstack:\n" << xerus::misc::get_call_stack());
~~~
The used macro will additionally include the source file and line as well as the function name in which the exception was thrown.

