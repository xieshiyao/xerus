# ------------------------------------------------------------------------------------------------------
#				Set the names of the resulting binary codes
# ------------------------------------------------------------------------------------------------------

# Names of the Library
LIB_NAME_SHARED = build/lib/libxerus.so
LIB_NAME_STATIC = build/lib/libxerus.a

# Name of the test executable
TEST_NAME = XerusTest

# xerus version from VERSION file
XERUS_VERSION = $(shell git describe --tags --always 2>/dev/null || cat VERSION)
DEBUG += -D XERUS_VERSION="$(XERUS_VERSION)"

XERUS_MAJOR_V = $(word 1, $(subst ., ,$(XERUS_VERSION)) )
ifneq (,$(findstring v, $(XERUS_MAJOR_V)))
	XERUS_MAJOR_V := $(strip $(subst v, ,$(XERUS_MAJOR_V)) )
endif
DEBUG += -DXERUS_VERSION_MAJOR=$(XERUS_MAJOR_V)
XERUS_MINOR_V = $(word 2, $(subst ., ,$(XERUS_VERSION)) )
DEBUG += -DXERUS_VERSION_MINOR=$(XERUS_MINOR_V)
XERUS_REVISION_V = $(word 3, $(subst ., ,$(XERUS_VERSION)) )
ifneq (,$(findstring -, $(XERUS_REVISION_V)))
	XERUS_COMMIT_V := $(word 2, $(subst -, ,$(XERUS_REVISION_V)) )
	XERUS_REVISION_V := $(word 1, $(subst -, ,$(XERUS_REVISION_V)) )
else
	XERUS_COMMIT_V = 0
endif
DEBUG += -DXERUS_VERSION_REVISION=$(XERUS_REVISION_V)
DEBUG += -DXERUS_VERSION_COMMIT=$(XERUS_COMMIT_V)


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

OTHER += -std=c++11			# Use the C++11 standard
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

# Small hack to get newlines...
define \n


endef

# Convinience variables
FLAGS = $(strip $(WARNINGS) $(OPTIMIZE) $(LOGGING) $(DEBUG) $(ADDITIONAL_INCLUDE) $(OTHER))
MINIMAL_DEPS = Makefile config.mk makeIncludes/general.mk makeIncludes/warnings.mk makeIncludes/optimization.mk



help:
	@printf "Possible make targets are:\n \
	\t\tall \t\t -- Build xerus both as a shared and a static library.\n \
	\t\tshared \t\t -- Build xerus as a shared library.\n \
	\t\tstatic \t\t -- Build xerus as a static library.\n \
	\t\ttest \t\t -- Build and run the xerus unit tests.\n \
	\t\tinstall \t -- Install the shared library and header files (may require root).\n \
	\t\t$(TEST_NAME) \t -- Only build the xerus unit tests.\n \
	\t\tclean \t\t -- Remove all object, library and executable files.\n"


opt:
	$(CXX) $(FLAGS) -Q --help=optimizers


warn:
	$(CXX) $(FLAGS) -Q --help=warnings


# Fake rule to create arbitary headers, to prevent errors if files are moved/renamed
%.h: 
	

all: test $(LIB_NAME_SHARED) $(LIB_NAME_STATIC)
	mkdir -p build/include/ 
	cp include/xerus.h build/include/
	cp -r include/xerus build/include/


shared: $(LIB_NAME_SHARED)


$(LIB_NAME_SHARED): $(MINIMAL_DEPS) $(LIB_SOURCES)
	mkdir -p $(dir $@)
	$(CXX) -shared -fPIC -Wl,-soname,libxerus.so $(FLAGS) -I include $(LIB_SOURCES) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -o $(LIB_NAME_SHARED) 


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


ifdef DESTDIR
	INSTALL_LIB_PATH = $(DESTDIR)
	INSTALL_HEADER_PATH = $(DESTDIR)/include
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
	


benchmark: $(MINIMAL_DEPS) $(LOCAL_HEADERS) benchmark.cxx $(LIB_NAME_STATIC)
	$(CXX) $(FLAGS) benchmark.cxx $(LIB_NAME_STATIC) $(SUITESPARSE) $(LAPACK_LIBRARIES) $(BLAS_LIBRARIES) $(CALLSTACK_LIBS) -lboost_filesystem -lboost_system -o Benchmark


# Build rule for normal lib objects
build/.libObjects/%.o: %.cpp $(MINIMAL_DEPS)
	mkdir -p $(dir $@) 
	$(CXX) -I include $< -c $(FLAGS) -MMD -o $@


# Build rule for test lib objects
build/.testObjects/%.o: %.cpp $(MINIMAL_DEPS)
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ -I include $< -c $(FLAGS) -MMD -o $@


# Build rule for benchmark objects
build/%.o: %.cpp $(MINIMAL_DEPS)
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ -I include $< -c $(FLAGS) -MMD -o $@


# Build rule for unit test objects
ifdef USE_GCC
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
build/.preCompileHeaders/xerus.h.gch: include/xerus.h $(MINIMAL_DEPS) .git/ORIG_HEAD
	mkdir -p $(dir $@)
	$(CXX) -D TEST_ $< -c $(FLAGS) -MMD -o $@


# dummy rule in case files were downloaded without git
.git/ORIG_HEAD:
	mkdir -p .git
	touch .git/ORIG_HEAD

