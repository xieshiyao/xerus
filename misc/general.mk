# Sets: CXX, LOCAL_HEADERS, LOCAL_HPP, TEST_HXX
# Uses: USE_CLANG, 

# Set Compiler used
ifdef USECLANG
	CXX = clang++
else
	CXX = g++
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

#Create Object list from general sources list and add path
# LIB_OBJECTS = $(LIB_SOURCES:%.cpp=.obj/%.o)
# TEST_OBJECTS = $(TEST_SOURCES:%.cpp=.obj/%.o)

