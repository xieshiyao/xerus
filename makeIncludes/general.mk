# Sets: CXX, CALLSTACK_LIBS, LOCAL_HEADERS, LOCAL_HPP, TEST_HXX
# Uses: USE_CLANG, 

# Set Compiler used
ifeq ($(strip $(CXX)),)
	ifdef USE_CLANG
		CXX = clang++
	else
		CXX = g++
	endif
endif


# include fancy_callstack specific libraries (binutils + dependencies)
ifdef NO_FANCY_CALLSTACK
	CALLSTACK_LIBS =
	DEBUG += -D NO_FANCY_CALLSTACK
else
	CALLSTACK_LIBS = -lbfd -liberty -lz -ldl 
endif


# Create list of local header files
LOCAL_HEADERS =  $(wildcard *.h)
LOCAL_HEADERS += $(wildcard */*.h)
LOCAL_HEADERS += $(wildcard */*/*.h)

# Create list of local .hpp files
LOCAL_HPP = $(wildcard *.hpp)
LOCAL_HPP += $(wildcard */*.hpp)
LOCAL_HPP += $(wildcard */*/*.hpp)

# Create list of local .hxx files
LOCAL_HXX = $(wildcard *.hxx)
LOCAL_HXX += $(wildcard */*.hxx)
LOCAL_HXX += $(wildcard */*/*.hxx)
