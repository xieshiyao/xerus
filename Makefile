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
LIB_SOURCES += $(wildcard src/xerus/*/*.cpp)

# Register all unit test source files 
UNIT_TEST_SOURCES = $(wildcard src/unitTests/*.cxx)

# Register all tutorial source files 
TUTORIAL_SOURCES = $(wildcard tutorials/*.cpp)

# Create lists of the corresponding objects and dependency files
LIB_OBJECTS = $(LIB_SOURCES:%.cpp=build/.libObjects/%.o)
LIB_DEPS    = $(LIB_SOURCES:%.cpp=build/.libObjects/%.d)

TEST_OBJECTS = $(LIB_SOURCES:%.cpp=build/.testObjects/%.o)
TEST_DEPS    = $(LIB_SOURCES:%.cpp=build/.testObjects/%.d)

UNIT_TEST_OBJECTS = $(UNIT_TEST_SOURCES:%.cxx=build/.unitTestObjects/%.o)
UNIT_TEST_DEPS    = $(UNIT_TEST_SOURCES:%.cxx=build/.unitTestObjects/%.d)

TUTORIALS 	= $(TUTORIAL_SOURCES:%.cpp=build/.tutorialObjects/%)
TUTORIAL_DEPS   = $(TUTORIAL_SOURCES:%.cpp=build/.tutorialObjects/%.d)

# Small hack to get newlines...
define \n


endef

# ------------------------------------------------------------------------------------------------------
#		Load the configurations provided by the user and set up general options
# ------------------------------------------------------------------------------------------------------
include config.mk

include makeIncludes/general.mk
include makeIncludes/warnings.mk
include makeIncludes/optimization.mk

# ------------------------------------------------------------------------------------------------------
#		  			Set additional compiler options                
# ------------------------------------------------------------------------------------------------------

# SUGGEST_ATTRIBUTES = TRUE 		# Tell the compiler to suggest attributes

# OPTIMIZE += -fprofile-generate	# Generate Profile output for optimization 
# OPTIMIZE += -fprofile-use		# Use a previously generated profile output

OTHER += -std=c++11			# Use old C++11 standard
OTHER += -D MISC_NAMESPACE=xerus	# All misc function shall live in xerus namespace

# ------------------------------------------------------------------------------------------------------
#					Load dependency files
# ------------------------------------------------------------------------------------------------------

-include $(LIB_DEPS)
-include $(TEST_DEPS)
-include $(UNIT_TEST_DEPS)
-include $(TUTORIAL_DEPS)
-include build/.preCompileHeaders/xerus.h.d



# ------------------------------------------------------------------------------------------------------
#					Make Rules      
# ------------------------------------------------------------------------------------------------------

# Convinience variables
FLAGS = $(strip $(WARNINGS) $(OPTIMIZE) $(LOGGING) $(DEBUG) $(ADDITIONAL_INCLUDE) $(OTHER))
MINIMAL_DEPS = Makefile makeIncludes/general.mk makeIncludes/warnings.mk makeIncludes/optimization.mk config.mk

help:
	@printf "Possible make targets are:\n \
	\t\tall \t\t\t -- Build xerus both as a shared and a static library.\n \
	\t\tshared \t\t -- Build xerus as a shared library.\n \
	\t\tstatic \t\t -- Build xerus as a static library.\n \
	\t\ttest \t\t\t -- Build and run the xerus unit tests.\n \
	\t\tinstall \t\t -- Install the shared library and header files (may require root).\n \
	\t\t$(TEST_NAME) \t\t -- Only build the xerus unit tests.\n \
	\t\tclean \t\t\t -- Remove all object, library and executable files.\n \
	\t\tbenchmark \t\t -- Build a loosly related benchmark program.\n \
	\t\tselectFunctions \t -- Performe tests to determine the best basic array functions to use on this machine.\n"

# Fake rule to create arbitary headers, to prevent errors if files are moved/renamed
%.h: 
	

all: test $(LIB_NAME_SHARED) $(LIB_NAME_STATIC)
	mkdir -p build/include/ 
	cp include/xerus.h build/include/
	cp -r include/xerus build/include/

shared: $(LIB_NAME_SHARED)

$(LIB_NAME_SHARED): $(MINIMAL_DEPS) $(LIB_SOURCES)
	mkdir -p $(dir $@)
	$(CXX) -shared -fPIC -Wl,-soname,libxerus.so $(FLAGS) -I include $(LIB_SOURCES) $(CALLSTACK_LIBS) -o $(LIB_NAME_SHARED) 

static: $(LIB_NAME_STATIC)

# Support non lto build for outdated systems
ifdef USE_LTO
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LIB_OBJECTS)
	mkdir -p $(dir $@)
	gcc-ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
else 
$(LIB_NAME_STATIC): $(MINIMAL_DEPS) $(LIB_OBJECTS)
	mkdir -p $(dir $@)
	ar rcs $(LIB_NAME_STATIC) $(LIB_OBJECTS)
endif

ifdef INSTALL_LIB_PATH
ifdef INSTALL_HEADER_PATH
install: $(LIB_NAME_SHARED)
	@printf "Installing libxerus.so to $(strip $(INSTALL_LIB_PATH)) and storing the header files in $(strip $(INSTALL_HEADER_PATH)).\n"
	mkdir -p $(INSTALL_LIB_PATH)
	mkdir -p $(INSTALL_HEADER_PATH)
	cp $(LIB_NAME_SHARED) $(INSTALL_LIB_PATH)
	cp include/xerus.h $(INSTALL_HEADER_PATH)
	cp -r include/xerus $(INSTALL_HEADER_PATH)
else
install:
	@printf "INSTALL_HEADER_PATH not set correctly. Cannot install xerus.\n"
endif
else
install:
	@printf "INSTALL_LIB_PATH not set correctly. Cannot install xerus.\n"
endif

$(TEST_NAME): $(MINIMAL_DEPS) $(UNIT_TEST_OBJECTS) $(TEST_OBJECTS)
	$(CXX) -D TEST_ $(FLAGS) $(UNIT_TEST_OBJECTS) $(TEST_OBJECTS) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o $(TEST_NAME)

test:  $(TEST_NAME)
	./$(TEST_NAME) all
	
fullTest: $(TUTORIALS) $(TEST_NAME)
	$(foreach x,$(TUTORIALS),./$(x)$(\n))
	./$(TEST_NAME) all

clean:
	rm -fr build
	-rm -f $(TEST_NAME)
	-rm -f include/xerus.h.gch
	

# selectFunctions: misc/preCompileSelector.cpp .obj/misc/stringUtilities.o .obj/misc/timeMeasure.o .obj/misc/namedLogger.o .obj/misc/blasLapackWrapper.o
# 	$(CXX) misc/preCompileSelector.cpp -std=c++11 -flto -fno-fat-lto-objects -flto-compression-level=0 --param ggc-min-heapsize=6442450 -Ofast -march=native \
# 	-D FULL_SELECTION_ libxerus.a $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o .obj/PreCompileSelector 
# 	.obj/PreCompileSelector

# 	
# benchmark: $(MINIMAL_DEPS) $(LOCAL_HEADERS) .obj/benchmark.o $(LIB_NAME_STATIC)
# 	$(CXX) -D CHECK_ $(FLAGS) .obj/benchmark.o $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) -o Benchmark
# 
# benchmarkTest: $(MINIMAL_DEPS) $(LOCAL_HEADERS) benchmark_tests.cxx $(LIB_NAME_STATIC)
# 	$(CXX) -D CHECK_ $(FLAGS) benchmark_tests.cxx $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o BenchmarkTest
# 

# Build rule for normal lib objects
build/.libObjects/%.o: %.cpp $(MINIMAL_DEPS)
	mkdir -p $(dir $@) 
	$(CXX) -I include $< -c $(FLAGS) -MMD -o $@

# Build rule for test lib objects
build/.testObjects/%.o: %.cpp $(MINIMAL_DEPS)
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ -I include $< -c $(FLAGS) -MMD -o $@

# Build rule for unit test objects
ifndef USE_CLANG
build/.unitTestObjects/%.o: %.cxx $(MINIMAL_DEPS) build/.preCompileHeaders/xerus.h.gch
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ -I build/.preCompileHeaders $< -c $(FLAGS) -MMD -o $@
else
build/.unitTestObjects/%.o: %.cxx $(MINIMAL_DEPS)
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ -I include $< -c $(FLAGS) -MMD -o $@
endif

# Build and execution rules for tutorials
build/.tutorialObjects/%: %.cpp $(MINIMAL_DEPS) $(LIB_NAME_STATIC)
	mkdir -p $(dir $@)
	$(CXX) -I include $< $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) $(FLAGS) -MMD -o $@

# Build rule for the preCompileHeader
build/.preCompileHeaders/xerus.h.gch: include/xerus.h $(MINIMAL_DEPS)
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ $< -c $(FLAGS) -MMD -o $@

