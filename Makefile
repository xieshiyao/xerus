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

# Register source files for the unit tests
TEST_SOURCES += xerusTest.cpp
TEST_SOURCES += unitTests/fullTensor_tests.cpp
TEST_SOURCES += misc/test.cpp

TEST_OBJECTS = $(TEST_SOURCES:%.cpp=.obj/%.o)


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
	$(CXX) -D TEST_ -shared -fPIC -Wl,-soname,libxerus.so $(FLAGS) $(LIB_SOURCES) -o $(LIB_NAME_SHARED) -liberty -lz -ldl -lbfd

# Support non lto build for outdated systems
ifdef LTO
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LIB_OBJECTS)
		gcc-ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
else 
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LIB_OBJECTS)     
		ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
endif

$(TEST_NAME): $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HXX) $(LIB_NAME_STATIC) $(TEST_OBJECTS) 
	$(CXX) -D TEST_ -D CHECK_ $(FLAGS) $(TEST_OBJECTS) $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -lbfd -liberty -lz -ldl -o $(TEST_NAME)

test:  $(TEST_NAME)
	./$(TEST_NAME) all
	
benchmark: $(MINIMAL_DEPS) $(LOCAL_HEADERS) .obj/benchmark.o $(LIB_NAME_STATIC)
	$(CXX) -D CHECK_ $(FLAGS) .obj/benchmark.o $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -o Benchmark

benchmarkTest: $(MINIMAL_DEPS) $(LOCAL_HEADERS) benchmark_tests.cxx $(LIB_NAME_STATIC)
	$(CXX) -D CHECK_ $(FLAGS) benchmark_tests.cxx $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -lbfd -liberty -lz -ldl -o BenchmarkTest

clean:
	-rm -f $(LIB_NAME) $(TEST_NAME) $(LIB_OBJECTS) $(TEST_OBJECTS)
	

selectFunctions: misc/preCompileSelector.cpp .obj/explodeString.o .obj/TimeMeasure.o .obj/timer.o .obj/log.o .obj/blasLapackWrapper.o
	$(CXX) includes/preCompileSelector.cpp -std=c++11 -I includes -I ../include -flto -fno-fat-lto-objects -flto-compression-level=0 --param ggc-min-heapsize=6442450 -Ofast -march=native \
	-fgcse-sm -fgcse-las -funswitch-loops -fipa-pta -fbranch-target-load-optimize -fsched-stalled-insns=100 -fsched-stalled-insns-dep=100 -fvariable-expansion-in-unroller \
	--param max-crossjump-edges=1000 --param max-grow-copy-bb-insns=10 --param max-delay-slot-insn-search=1000 --param max-delay-slot-live-search=1000 --param max-gcse-memory=6442450 \
	--param max-pending-list-length=1000 --param max-modulo-backtrack-attempts=1000 --param max-inline-insns-single=4000 --param max-inline-insns-auto=400 --param large-function-insns=10000 \
	--param large-function-growth=200 --param inline-unit-growth=100 --param ipcp-unit-growth=100 --param max-reload-search-insns=500 --param max-cselib-memory-locations=5000 \
	--param max-sched-ready-insns=1000 --param max-sched-region-blocks=100 --param max-pipeline-region-blocks=150 --param max-sched-region-insns=1000 --param max-pipeline-region-insns=2000 \
	--param selsched-max-lookahead=500 --param max-partial-antic-length=0 --param sccvn-max-scc-size=100000 --param sccvn-max-alias-queries-per-access=10000 --param ira-max-loops-num=1000 \
	--param ira-max-conflict-table-size=20000 --param loop-invariant-max-bbs-in-loop=100000 --param loop-max-datarefs-for-datadeps=10000 --param max-vartrack-size=0 --param max-vartrack-expr-depth=120 \
	-freciprocal-math -fmerge-all-constants -D FULL_SELECTION_ .obj/explodeString.o .obj/TimeMeasure.o .obj/timer.o .obj/log.o .obj/blasLapackWrapper.o $(EXTRA_LAPACK) -o includes/.obj/PreCompileSelector
	includes/.obj/PreCompileSelector
	

#Compile local source files - depend on all Headers and local directorys
.obj/%.o: %.cpp $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HPP)
	mkdir -p $(dir $@)
	$(CXX) $< -c -D TEST_ -D CHECK_ $(FLAGS) -o $@

#Compile local source files - depend on all Headers and local directorys
.obj/%.o: %.cxx $(MINIMAL_DEPS) $(LOCAL_HEADERS) $(LOCAL_HPP)
	mkdir -p $(dir $@)
	$(CXX) $< -c -D TEST_ -D CHECK_ $(FLAGS) -o $@
