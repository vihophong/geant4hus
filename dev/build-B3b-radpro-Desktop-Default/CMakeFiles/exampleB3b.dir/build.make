# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lubuntu/phong/geant4hus/dev/B3b-radpro

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default

# Include any dependencies generated for this target.
include CMakeFiles/exampleB3b.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/exampleB3b.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/exampleB3b.dir/flags.make

CMakeFiles/exampleB3b.dir/exampleB3b.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/exampleB3b.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/exampleB3b.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/exampleB3b.dir/exampleB3b.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/exampleB3b.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/exampleB3b.cc

CMakeFiles/exampleB3b.dir/exampleB3b.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/exampleB3b.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/exampleB3b.cc > CMakeFiles/exampleB3b.dir/exampleB3b.cc.i

CMakeFiles/exampleB3b.dir/exampleB3b.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/exampleB3b.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/exampleB3b.cc -o CMakeFiles/exampleB3b.dir/exampleB3b.cc.s

CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.requires

CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.provides: CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.provides

CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.provides.build: CMakeFiles/exampleB3b.dir/exampleB3b.cc.o


CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3DetectorConstruction.cc

CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3DetectorConstruction.cc > CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.i

CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3DetectorConstruction.cc -o CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.s

CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o


CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PhysicsList.cc

CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PhysicsList.cc > CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.i

CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PhysicsList.cc -o CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.s

CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o


CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PrimaryGeneratorAction.cc

CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PrimaryGeneratorAction.cc > CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.i

CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3PrimaryGeneratorAction.cc -o CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.s

CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o


CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3StackingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3StackingAction.cc

CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3StackingAction.cc > CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.i

CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3StackingAction.cc -o CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.s

CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o


CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bActionInitialization.cc

CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bActionInitialization.cc > CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.i

CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bActionInitialization.cc -o CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.s

CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o


CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRun.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRun.cc

CMakeFiles/exampleB3b.dir/src/B3bRun.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3bRun.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRun.cc > CMakeFiles/exampleB3b.dir/src/B3bRun.cc.i

CMakeFiles/exampleB3b.dir/src/B3bRun.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3bRun.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRun.cc -o CMakeFiles/exampleB3b.dir/src/B3bRun.cc.s

CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o


CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o: CMakeFiles/exampleB3b.dir/flags.make
CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o: /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o -c /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRunAction.cc

CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRunAction.cc > CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.i

CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lubuntu/phong/geant4hus/dev/B3b-radpro/src/B3bRunAction.cc -o CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.s

CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.requires

CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.provides: CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3b.dir/build.make CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.provides

CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.provides.build: CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o


# Object files for target exampleB3b
exampleB3b_OBJECTS = \
"CMakeFiles/exampleB3b.dir/exampleB3b.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o" \
"CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o"

# External object files for target exampleB3b
exampleB3b_EXTERNAL_OBJECTS =

exampleB3b: CMakeFiles/exampleB3b.dir/exampleB3b.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o
exampleB3b: CMakeFiles/exampleB3b.dir/build.make
exampleB3b: /home/lubuntu/geant4install/lib/libG4Tree.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4GMocren.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4visHepRep.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4RayTracer.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4VRML.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4OpenGL.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4gl2ps.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4interfaces.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4persistency.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4error_propagation.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4readout.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4physicslists.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4parmodels.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4FR.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4vis_management.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4modeling.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libXmu.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libXext.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libXt.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libSM.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libICE.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libX11.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libXm.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libGLU.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libGL.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libQtGui.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libQtCore.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libxerces-c.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4run.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4event.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4tracking.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4processes.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4analysis.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4zlib.so
exampleB3b: /usr/lib/x86_64-linux-gnu/libexpat.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4digits_hits.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4track.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4particles.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4geometry.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4materials.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4graphics_reps.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4intercoms.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4global.so
exampleB3b: /home/lubuntu/geant4install/lib/libG4clhep.so
exampleB3b: CMakeFiles/exampleB3b.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable exampleB3b"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/exampleB3b.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/exampleB3b.dir/build: exampleB3b

.PHONY : CMakeFiles/exampleB3b.dir/build

CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/exampleB3b.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3DetectorConstruction.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3PhysicsList.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3PrimaryGeneratorAction.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3StackingAction.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3bActionInitialization.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3bRun.cc.o.requires
CMakeFiles/exampleB3b.dir/requires: CMakeFiles/exampleB3b.dir/src/B3bRunAction.cc.o.requires

.PHONY : CMakeFiles/exampleB3b.dir/requires

CMakeFiles/exampleB3b.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/exampleB3b.dir/cmake_clean.cmake
.PHONY : CMakeFiles/exampleB3b.dir/clean

CMakeFiles/exampleB3b.dir/depend:
	cd /home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lubuntu/phong/geant4hus/dev/B3b-radpro /home/lubuntu/phong/geant4hus/dev/B3b-radpro /home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default /home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default /home/lubuntu/phong/geant4hus/dev/build-B3b-radpro-Desktop-Default/CMakeFiles/exampleB3b.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/exampleB3b.dir/depend

