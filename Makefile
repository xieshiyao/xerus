# ------------------------------------------------------------------------------------------------------
#				Set the names of the resulting binary codes
# ------------------------------------------------------------------------------------------------------

# Names of the Library
LIB_NAME_SHARED = libxerus.so
LIB_NAME_STATIC = libxerus.a

# Name of the test executable
TEST_NAME = XerusTest


# ------------------------------------------------------------------------------------------------------
#				Register source files for the xerus library      
# ------------------------------------------------------------------------------------------------------

# Register source files for the library
LIB_SOURCES += xerus.cpp
LIB_SOURCES += misc/timeMeasure.cpp
LIB_SOURCES += misc/stringUtilities.cpp
LIB_SOURCES += misc/namedLogger.cpp
LIB_SOURCES += misc/blasLapackWrapper.cpp
LIB_SOURCES += misc/callStack.cpp
LIB_SOURCES += misc/simpleNumerics.cpp

LIB_OBJECTS = $(LIB_SOURCES:%.cpp=.obj/%.o)
TEST_LIB_OBJECTS = $(LIB_SOURCES:%.cpp=.testObj/%.o)

# Register source files for the unit tests
TEST_SOURCES += xerusTest.cpp
TEST_SOURCES += unitTests/fullTensor_tests.cpp
TEST_SOURCES += misc/test.cpp

TEST_OBJECTS = $(TEST_SOURCES:%.cpp=.testObj/%.o)


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
include misc/general.mk
include misc/warnings.mk
include misc/optimization.mk

FLAGS = $(strip $(LOGGING) $(DEBUG) $(WARNINGS) $(OPTIMIZE) $(ADDITIONAL_INCLUDE) $(OTHER))
MINIMAL_DEPS = Makefile misc/general.mk misc/warnings.mk misc/optimization.mk config.mk
# ------------------------------------------------------------------------------------------------------
#					Make Rules      
# ------------------------------------------------------------------------------------------------------
help:
	@printf "Possible make targets are:\n \
	\t\t$(LIB_NAME_SHARED) \t\t -- Build Xerus as a shared library.\n \
	\t\t$(LIB_NAME_STATIC) \t\t -- Build Xerus as a static library.\n \
	\t\t$(TEST_NAME) \t\t -- Build the Xerus unit tests.\n \
	\t\ttest \t\t\t -- Run the Xerus unit tests.\n \
	\t\tbenchmark \t\t -- Build a loosly related benchmark program.\n \
	\t\tclean \t\t\t -- Remove all object and library and executable files.\n \
	\t\tselectFunctions \t -- Performe tests to determine the best basic array functions to use on this machine.\n"
	
debug:
	@printf "Some variables we have set:\n \
	\t\t LIB_NAME_SHARED \t == \t $(LIB_NAME_SHARED)\n \
	\t\t LIB_NAME_STATIC \t == \t $(LIB_NAME_STATIC)\n \
	\t\t TEST_NAME \t == \t $(TEST_NAME)\n \
	\t\t LIB_SOURCES \t == \t $(LIB_SOURCES)\n \
	\t\t LIB_OBJECTS \t == \t $(LIB_OBJECTS)\n \
	\t\t TEST_SOURCES \t == \t $(TEST_SOURCES)\n \
	\t\t TEST_OBJECTS \t == \t $(TEST_OBJECTS)\n \
	\t\t LOCAL_HEADERS \t == \t $(LOCAL_HEADERS)\n \
	\t\t LOCAL_HPP \t == \t $(LOCAL_HPP)\n \
	\t\t LOCAL_HXX \t == \t $(LOCAL_HXX)\n \
	\t\t MINIMAL_DEPS \t == \t $(MINIMAL_DEPS)\n"

$(LIB_NAME_SHARED): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LIB_SOURCES)
	$(CXX) -shared -fPIC -Wl,-soname,libxerus.so $(FLAGS) $(LIB_SOURCES) -o $(LIB_NAME_SHARED) -liberty -lz -ldl -lbfd

# Support non lto build for outdated systems
ifdef USE_LTO
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LIB_OBJECTS)
	gcc-ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
else 
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LIB_OBJECTS)     
	ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
endif

$(TEST_NAME): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HXX) $(TEST_LIB_OBJECTS) $(TEST_OBJECTS) 
	$(CXX) -D TEST_ $(FLAGS) $(TEST_LIB_OBJECTS) $(TEST_OBJECTS) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -lbfd -liberty -lz -ldl -o $(TEST_NAME)

test:  $(TEST_NAME)
	./$(TEST_NAME) all
	
benchmark: $(MINIMAL_DEPS) $(LOCAL_HEADERS) .obj/benchmark.o $(LIB_NAME_STATIC)
	$(CXX) -D CHECK_ $(FLAGS) .obj/benchmark.o $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -o Benchmark

benchmarkTest: $(MINIMAL_DEPS) $(LOCAL_HEADERS) benchmark_tests.cxx $(LIB_NAME_STATIC)
	$(CXX) -D CHECK_ $(FLAGS) benchmark_tests.cxx $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -lbfd -liberty -lz -ldl -o BenchmarkTest

clean:
	-rm -f $(LIB_NAME_STATIC) $(LIB_NAME_SHARED) $(TEST_NAME) $(LIB_OBJECTS) $(TEST_LIB_OBJECTS) $(TEST_OBJECTS) .obj/PreCompileSelector
	

selectFunctions: misc/preCompileSelector.cpp .obj/misc/stringUtilities.o .obj/misc/timeMeasure.o .obj/misc/namedLogger.o .obj/misc/blasLapackWrapper.o
	$(CXX) misc/preCompileSelector.cpp -std=c++11 -flto -fno-fat-lto-objects -flto-compression-level=0 --param ggc-min-heapsize=6442450 -Ofast -march=native \
	-D FULL_SELECTION_ libxerus.a $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -lbfd -liberty -lz -ldl -o .obj/PreCompileSelector 
	.obj/PreCompileSelector
	
# Compile sources to test object files
.testObj/%.o: %.cpp $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HPP)
	mkdir -p $(dir $@)
	$(CXX) $< -c -D TEST_ $(FLAGS) -o $@

.testObj/%.o: %.cxx $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HPP) $(LOCAL_HXX)
	mkdir -p $(dir $@)
	$(CXX) $< -c -D TEST_ $(FLAGS) -o $@

#Compile local source files - depend on all Headers and local directorys
.obj/%.o: %.cpp $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HPP)
	mkdir -p $(dir $@)
	$(CXX) $< -c $(FLAGS) -o $@

