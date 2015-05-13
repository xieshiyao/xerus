# ------------------------------------------------------------------------------------------------------
#				Set the names of the resulting binary codes
# ------------------------------------------------------------------------------------------------------

# Names of the Library
LIB_NAME_SHARED = build/lib/libxerus.so
LIB_NAME_STATIC = build/lib/libxerus.a

# Name of the test executable
TEST_NAME = XerusTest


# ------------------------------------------------------------------------------------------------------
#				Register source files for the xerus library      
# ------------------------------------------------------------------------------------------------------

# Register all library source files 
LIB_SOURCES = $(wildcard src/xerus/*.cpp)
LIB_SOURCES += $(wildcard src/xerus/algorithm/*.cpp)
LIB_SOURCES_NO_DIR = $(notdir $(LIB_SOURCES))

# Register all misc source files 
MISC_SOURCES = $(wildcard src/misc/*.cpp)
MISC_SOURCES_NO_DIR = $(notdir $(MISC_SOURCES))

# Register all unit test source files 
UNIT_TEST_SOURCES = $(wildcard src/unitTests/*.cpp)
UNIT_TEST_SOURCES_NO_DIR = $(notdir $(UNIT_TEST_SOURCES))

LIB_OBJECTS = $(LIB_SOURCES_NO_DIR:%.cpp=build/lib/.libObjects/%.o)
LIB_DEPS    = $(LIB_SOURCES_NO_DIR:%.cpp=build/lib/.libObjects/%.d)

MISC_OBJECTS = $(MISC_SOURCES_NO_DIR:%.cpp=build/lib/.miscObjects/%.o)
MISC_DEPS    = $(MISC_SOURCES_NO_DIR:%.cpp=build/lib/.miscObjects/%.d)

TEST_LIB_OBJECTS = $(LIB_SOURCES_NO_DIR:%.cpp=build/lib/.testLibObjects/%.o)
TEST_LIB_DEPS    = $(LIB_SOURCES_NO_DIR:%.cpp=build/lib/.testLibObjects/%.d)

TEST_MISC_OBJECTS = $(MISC_SOURCES_NO_DIR:%.cpp=build/lib/.testMiscObjects/%.o)
TEST_MISC_DEPS    = $(MISC_SOURCES_NO_DIR:%.cpp=build/lib/.testMiscObjects/%.d)

UNIT_TEST_OBJECTS = $(UNIT_TEST_SOURCES_NO_DIR:%.cpp=build/lib/.unitTestObjects/%.o)
UNIT_TEST_DEPS    = $(UNIT_TEST_SOURCES_NO_DIR:%.cpp=build/lib/.unitTestObjects/%.d)


# ------------------------------------------------------------------------------------------------------
#				Load the configurations provided by the user  
# ------------------------------------------------------------------------------------------------------
include config.mk


# ------------------------------------------------------------------------------------------------------
#		  			Setup compiler options                
# ------------------------------------------------------------------------------------------------------

# SUGGEST_ATTRIBUTES = TRUE 		# Tell the compiler to suggest attributes

  
# OPTIMIZE += -fprofile-generate	# Generate Profile output for optimization 
# OPTIMIZE += -fprofile-use		# Use a previously generated profile output

OTHER += -I misc 			# Add search Path for header files
OTHER += -std=c++11			# Use old C++11 standard
OTHER += -D MISC_NAMESPACE=xerus	# All misc function shall live in xerus namespace


# ------------------------------------------------------------------------------------------------------
#					Setup general compiler options      
# ------------------------------------------------------------------------------------------------------
include makeIncludes/general.mk
include makeIncludes/warnings.mk
include makeIncludes/optimization.mk

FLAGS = $(strip $(LOGGING) $(DEBUG) $(WARNINGS) $(OPTIMIZE) $(ADDITIONAL_INCLUDE) $(OTHER))
MINIMAL_DEPS = Makefile makeIncludes/general.mk makeIncludes/warnings.mk makeIncludes/optimization.mk config.mk build/.structure
# ------------------------------------------------------------------------------------------------------
#					Load dependency files
# ------------------------------------------------------------------------------------------------------

-include $(LIB_DEPS)
-include $(MISC_DEPS)
-include $(TEST_LIB_DEPS)
-include $(TEST_MISC_DEPS)
-include $(UNIT_TEST_DEPS)

# ------------------------------------------------------------------------------------------------------
#					Make Rules      
# ------------------------------------------------------------------------------------------------------
help:
	@printf "Possible make targets are:\n \
	\t\tall \t\t\t -- Build xerus both as a shared and a static library.\n \
	\t\t$(LIB_NAME_SHARED) \t\t -- Build xerus as a shared library.\n \
	\t\t$(LIB_NAME_STATIC) \t\t -- Build xerus as a static library.\n \
	\t\ttest \t\t\t -- Build and run the xerus unit tests.\n \
	\t\tinstall \t\t -- Install the shared library and header files (may require root).\n \
	\t\t$(TEST_NAME) \t\t -- Only build the xerus unit tests.\n \
	\t\tbenchmark \t\t -- Build a loosly related benchmark program.\n \
	\t\tclean \t\t\t -- Remove all object and library and executable files.\n \
	\t\tselectFunctions \t -- Performe tests to determine the best basic array functions to use on this machine.\n \
	\t\tdebug \t\t\t -- Print information on the Makefile variables set.\n"
	
debug:
	@printf "Some variables we have set:\n \
	\t\t LIB_NAME_SHARED \t == \t $(LIB_NAME_SHARED)\n\n \
	\t\t LIB_NAME_STATIC \t == \t $(LIB_NAME_STATIC)\n\n \
	\t\t TEST_NAME \t\t == \t $(TEST_NAME)\n\n \
	\t\t LIB_SOURCES \t\t == \t $(LIB_SOURCES)\n\n \
	\t\t LIB_OBJECTS \t\t == \t $(LIB_OBJECTS)\n\n \
	\t\t LIB_DEPS \t\t == \t $(LIB_DEPS)\n\n \
	\t\t TEST_OBJECTS \t\t == \t $(TEST_OBJECTS)\n\n \
	\t\t LOCAL_HEADERS \t\t == \t $(LOCAL_HEADERS)\n\n \
	\t\t LOCAL_HPP \t\t == \t $(LOCAL_HPP)\n\n \
	\t\t LOCAL_HXX \t\t == \t $(LOCAL_HXX)\n\n \
	\t\t MINIMAL_DEPS \t\t == \t $(MINIMAL_DEPS)\n\n"

build/.structure:
	mkdir -p build build/include build/include build/include/.preCompiledHeaders build/lib build/lib/.libObjects build/lib/.miscObjects build/lib/.testLibObjects build/lib/.testMiscObjects build/lib/.unitTestObjects
	touch build/.structure

all: $(LIB_NAME_SHARED) $(LIB_NAME_STATIC)
	cp include/xerus.h build/include/
	cp -r include/xerus build/include/
	cp -r include/misc build/include/

$(LIB_NAME_SHARED): $(MINIMAL_DEPS) $(LIB_SOURCES) $(MISC_SOURCES)
	$(CXX) -shared -fPIC -Wl,-soname,libxerus.so $(FLAGS) $(LIB_SOURCES) $(MISC_SOURCES) $(CALLSTACK_LIBS) -o $(LIB_NAME_SHARED) 

# Support non lto build for outdated systems
ifdef USE_LTO
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LIB_OBJECTS) $(MISC_OBJECTS)
	gcc-ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS) $(MISC_OBJECTS)
else 
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LIB_OBJECTS) $(MISC_OBJECTS) 
	ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS) $(MISC_OBJECTS)
endif

install:
	@printf "Sorry not yet supported\n" # TODO

$(TEST_NAME): $(MINIMAL_DEPS) $(UNIT_TEST_OBJECTS) $(TEST_LIB_OBJECTS) $(TEST_MISC_OBJECTS)
	$(CXX) -D TEST_ $(FLAGS) $(UNIT_TEST_OBJECTS) $(TEST_LIB_OBJECTS) $(TEST_MISC_OBJECTS) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o $(TEST_NAME)

test:  $(TEST_NAME)
	./$(TEST_NAME) all
# 	
# benchmark: $(MINIMAL_DEPS) $(LOCAL_HEADERS) .obj/benchmark.o $(LIB_NAME_STATIC)
# 	$(CXX) -D CHECK_ $(FLAGS) .obj/benchmark.o $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -o Benchmark
# 
# benchmarkTest: $(MINIMAL_DEPS) $(LOCAL_HEADERS) benchmark_tests.cxx $(LIB_NAME_STATIC)
# 	$(CXX) -D CHECK_ $(FLAGS) benchmark_tests.cxx $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o BenchmarkTest
# 
clean:
	rm -fr build
	-rm -f $(TEST_NAME)
	-rm -f include/xerus.h.gch
	
# 
# selectFunctions: misc/preCompileSelector.cpp .obj/misc/stringUtilities.o .obj/misc/timeMeasure.o .obj/misc/namedLogger.o .obj/misc/blasLapackWrapper.o
# 	$(CXX) misc/preCompileSelector.cpp -std=c++11 -flto -fno-fat-lto-objects -flto-compression-level=0 --param ggc-min-heapsize=6442450 -Ofast -march=native \
# 	-D FULL_SELECTION_ libxerus.a $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o .obj/PreCompileSelector 
# 	.obj/PreCompileSelector

build/lib/.libObjects%.o: src/xerus%.cpp $(MINIMAL_DEPS)
	$(CXX) $< -c $(FLAGS) -MMD -o $@

build/lib/.miscObjects%.o: src/misc%.cpp $(MINIMAL_DEPS)
	$(CXX) $< -c $(FLAGS) -MMD -o $@

build/lib/.testLibObjects%.o: src/xerus%.cpp $(MINIMAL_DEPS)
	$(CXX) -D TEST_ $< -c $(FLAGS) -MMD -o $@

build/lib/.testMiscObjects%.o: src/misc%.cpp $(MINIMAL_DEPS)
	$(CXX) -D TEST_ $< -c $(FLAGS) -MMD -o $@

build/lib/.unitTestObjects%.o: src/unitTests%.cpp $(MINIMAL_DEPS) include/xerus.h.gch
	$(CXX) -D TEST_ $< -c $(FLAGS) -MMD -o $@
	
include/xerus.h.gch: include/xerus.h $(MINIMAL_DEPS)
	$(CXX) -D TEST_ $< $(FLAGS) -MMD -o $@

