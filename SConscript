##
# @file SConscript
# @brief Provides functionality for building the executable.
Import('env')

boostLibs = [
             'boost_filesystem',
             'boost_program_options',
             'boost_system',
             ]

srcFiles = [
            'ProgramOptions.cpp',
            'driver.cpp',
            ]

built = False
# Do not modify the provided build environment
buildEnv = env.Clone()
# Check if all the library files can be located
conf = Configure(buildEnv)
for lib in boostLibs:
  if not conf.CheckLib(lib, language='C++'):
    Return('built')
buildEnv = conf.Finish()

buildEnv.Program(target=buildEnv.targetName, source=srcFiles)

buildEnv.Install(buildEnv.topDir, buildEnv.targetName)
built = True
Return('built')
