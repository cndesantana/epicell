###########################################################
# Configuration file for the compilation of Epicell code,
# using the SConstruct library.
# IT IS NOT RECOMMENDED TO MODIFY THIS FILE.
# Compilation should be personalized by adjusting the 
# Makefile in the directory of the main source files.
# See Palabos examples for sample Makefiles.
###########################################################

import os
import sys
import glob

argdict = dict(ARGLIST)

# Read input parameters
epicellRoot   = argdict['epicellRoot']
projectFiles  = Split(argdict['projectFiles'])
optimize      = argdict['optimize'].lower() == 'true'
debug         = argdict['debug'].lower() == 'true'
profile       = argdict['profile'].lower() == 'true'
MPIparallel   = argdict['MPIparallel'].lower() == 'true'
serialCXX     = argdict['serialCXX']
parallelCXX   = argdict['parallelCXX']
compileFlags  = Split(argdict['compileFlags'])
linkFlags     = Split(argdict['linkFlags'])
optimFlags    = Split(argdict['optimFlags'])
debugFlags    = Split(argdict['debugFlags'])
profileFlags  = Split(argdict['profileFlags'])
libraryPaths  = Split(argdict['libraryPaths'])
includePaths  = Split(argdict['includePaths'])
libraries     = Split(argdict['libraries'])

# Read the optional input parameters
try:
    dynamicLibrary = argdict['dynamicLibrary'].lower() == 'true'
except:
    dynamicLibrary = False

try:
    srcPaths = Split(argdict['srcPaths'])
except:
    srcPaths = []

flags = compileFlags
allPaths = [epicellRoot+'/src'] + [epicellRoot+'/externalLibraries'] + includePaths

if optimize:
    flags.append(optimFlags)

if debug:
    flags.append(debugFlags)
    flags.append('-DEPC_DEBUG')

if profile:
    flags.append(profileFlags)
    linkFlags.append(profileFlags)

if MPIparallel:
    compiler = parallelCXX
    flags.append('-DEPC_MPI_PARALLEL')
else:
    compiler = serialCXX

env = Environment ( ENV       = os.environ,
                    CXX       = compiler,
                    CXXFLAGS  = flags,
                    LINKFLAGS = linkFlags,
                    CPPPATH   = allPaths
                  )

if dynamicLibrary:
    LibraryGen = env.SharedLibrary
else:
    LibraryGen = env.Library


sourceFiles = []
for srcDir in glob.glob(epicellRoot+'/src/*'):
    sourceFiles.extend(glob.glob(srcDir+'/*.cpp'))

for srcDir in srcPaths:
    sourceFiles.extend(glob.glob(srcDir+'/*.cpp'))

sourceFiles.extend(glob.glob(epicellRoot+'/externalLibraries/tinyxml/*.cpp'));

if MPIparallel:
    epicell_library = LibraryGen( target  = epicellRoot+'/lib/plb_mpi',
                                  source  = sourceFiles )
else:
    epicell_library = LibraryGen( target  = epicellRoot+'/lib/plb',
                                  source  = sourceFiles )

local_objects = env.Object(source = projectFiles)

all_objects = local_objects + epicell_library

env.Program(all_objects, LIBS=libraries, LIBPATH=libraryPaths)
