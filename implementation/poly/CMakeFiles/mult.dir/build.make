# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/cusgadmin/Desktop/poly

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/cusgadmin/Desktop/poly

# Include any dependencies generated for this target.
include CMakeFiles/mult.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mult.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mult.dir/flags.make

CMakeFiles/mult.dir/multivariate_test.cpp.o: CMakeFiles/mult.dir/flags.make
CMakeFiles/mult.dir/multivariate_test.cpp.o: multivariate_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cusgadmin/Desktop/poly/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mult.dir/multivariate_test.cpp.o"
	/Applications/Xcode9.4.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mult.dir/multivariate_test.cpp.o -c /Users/cusgadmin/Desktop/poly/multivariate_test.cpp

CMakeFiles/mult.dir/multivariate_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mult.dir/multivariate_test.cpp.i"
	/Applications/Xcode9.4.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/cusgadmin/Desktop/poly/multivariate_test.cpp > CMakeFiles/mult.dir/multivariate_test.cpp.i

CMakeFiles/mult.dir/multivariate_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mult.dir/multivariate_test.cpp.s"
	/Applications/Xcode9.4.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/cusgadmin/Desktop/poly/multivariate_test.cpp -o CMakeFiles/mult.dir/multivariate_test.cpp.s

# Object files for target mult
mult_OBJECTS = \
"CMakeFiles/mult.dir/multivariate_test.cpp.o"

# External object files for target mult
mult_EXTERNAL_OBJECTS =

mult: CMakeFiles/mult.dir/multivariate_test.cpp.o
mult: CMakeFiles/mult.dir/build.make
mult: CMakeFiles/mult.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/cusgadmin/Desktop/poly/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mult"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mult.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mult.dir/build: mult

.PHONY : CMakeFiles/mult.dir/build

CMakeFiles/mult.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mult.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mult.dir/clean

CMakeFiles/mult.dir/depend:
	cd /Users/cusgadmin/Desktop/poly && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/cusgadmin/Desktop/poly /Users/cusgadmin/Desktop/poly /Users/cusgadmin/Desktop/poly /Users/cusgadmin/Desktop/poly /Users/cusgadmin/Desktop/poly/CMakeFiles/mult.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mult.dir/depend

