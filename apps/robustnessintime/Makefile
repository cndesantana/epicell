##########################################################################
## Makefile.
##
## The present Makefile is a pure configuration file, in which 
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the Palabos source code
#epicellRoot  = $(EPICELL_ROOT)
epicellRoot  = /home/cdesantana/Data/GitHubRepo/epicell

projectFiles = robustnessintime.cpp

# Set optimization flags on/off
optimize     = false
# Set debug mode and debug flags on/off
debug        = true
# Set profiling flags on/off
profile      = false
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = false

# Path to external libraries (other than Palabos)
libraryPaths =
# Path to inlude directories (other than Palabos)
includePaths =
# Dynamic and static libraries (other than Palabos)
libraries    =

# Compiler to use without MPI parallelism
#serialCXX    = clang++
serialCXX    = g++
# Compiler to use with MPI parallelism
parallelCXX  = mpicxx
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
compileFlags = -Wall -std=c++0x 
# General linker flags (don't put library includes into this flag)
linkFlags    =
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
# Compiler flags to use when debug mode is on
debugFlags   = -g3
# Compiler flags to use when profile mode is on
profileFlags = -pg


##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(epicellRoot)/scons/scons.py -j 6 -f $(epicellRoot)/SConstruct

SConsArgs = epicellRoot=$(epicellRoot) \
            projectFiles="$(projectFiles)" \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
            profileFlags="$(profileFlags)" \
            libraryPaths="$(libraryPaths)" \
            includePaths="$(includePaths)" \
            libraries="$(libraries)"

compile:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(epicellRoot) -name '*~'`
