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
# Check if all the library files can be located
conf = Configure(env)
for lib in boostLibs:
  if not conf.CheckLib(lib, language='C++'):
    Return('built')
env = conf.Finish()

env.Program(target=env.targetName, source=srcFiles)

env.Install(env.topDir, env.targetName)
built = True
Return('built')
