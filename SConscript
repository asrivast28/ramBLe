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

allLibs = env.get('LIBS', [])

env.Program(target=env.targetName, source=srcFiles, LIBS=allLibs+boostLibs)

env.Install(env.topDir, env.targetName)
