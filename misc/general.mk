# Sets: CXX, LOCAL_HEADERS, LOCAL_HPP, TEST_HXX
# Uses: USE_CLANG, 

# Set Compiler used
ifdef USE_CLANG
	CXX = clang++
else
	CXX = g++
endif


# include fancy_callstack specific libraries (binutils + dependencies)
ifdef FANCY_CALLSTACK
	CALLSTACK_LIBS = -lbfd -liberty -lz -ldl
	OTHER += -D FANCY_CALLSTACK_
else
	CALLSTACK_LIBS = 
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
