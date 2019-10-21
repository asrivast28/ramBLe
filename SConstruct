##
# @file SConstruct
# @brief Provides general functionality for building.
import os
import platform

topDir = os.path.abspath(os.getcwd())
cpp = None
buildDir = None
targetName = 'csl'
testName = 'testcsl'


cppPaths = [
            topDir,
            ]

cppFlags = [
            ]

cppDefs = [
           ]

libPaths = [
            ]

allLibs = [
            ]

linkFlags = [
             ]

# Flag for building in debug mode. Defaults to release build.
releaseBuild = ARGUMENTS.get('DEBUG', 0) in [0, '0']
# Comma separated list of non-default library locations
localLibPaths = ARGUMENTS.get('LOCALLIBS')
if localLibPaths is not None:
  libPaths.extend(localLibPaths.split())
# Comma separated list of non-default include directory locations
localIncludePaths = ARGUMENTS.get('LOCALINCLUDES')
if localIncludePaths is not None:
  cppPaths.extend(localIncludePaths.split())
# Location of the SABNAtk directory.
sabnatkDir = ARGUMENTS.get('SABNATK', os.path.join(topDir, 'SABNAtk'))
if os.path.exists(sabnatkDir):
  cppPaths.append(os.path.join(sabnatkDir, 'include'))


if platform.system() in ['Darwin', 'Linux']:
  cpp = 'g++'
  cppFlags.extend([
              '-Wall',
              '-std=c++14',
              ])
  if releaseBuild:
    cppFlags.append('-O3')
    cppDefs.append('NDEBUG')
  else:
    cppFlags.append('-g')

# For OS X
if platform.system() == 'Darwin':
  libPaths.append('/opt/local/lib')
  cppPaths.append('/opt/local/include')

if releaseBuild:
  buildDir = 'release'
else:
  buildDir = 'debug'
  targetName += '_debug'
  testName += '_debug'

enableLogging = ARGUMENTS.get('LOGGING', None)
if enableLogging is not None:
  # use the user provided option for enabling/disabling logging
  enableLogging = False if enableLogging in [0, '0'] else True
else:
  # otherwise, enable logging only in debug build
  enableLogging = not releaseBuild

if enableLogging:
  # logging libraries
  cppDefs.append('LOGGING')
  allLibs.extend(['boost_log', 'boost_system', 'pthread'])

# Flag for enabling profiling
enableProfile = ARGUMENTS.get('PROFILE', 0) not in [0, '0']
if enableProfile:
  if platform.system() == 'Linux':
    cppFlags.append('-pg')
    linkFlags.append('-pg')
    targetName += '_profile'
    testName += '_profile'
  else:
    print("WARNING: Profiling is not supported on", platform.system())

env = Environment(ENV=os.environ, CXX=cpp, CXXFLAGS=cppFlags, CPPPATH=cppPaths, CPPDEFINES=cppDefs, LIBPATH=libPaths, LINKFLAGS=linkFlags)
conf = Configure(env)
# Check if the initial build environment works
if not conf.CheckCXX():
  Exit(1)
# Check for SABNAtk specific functions and build options
sabnatkBuiltins = ['clzll', 'ffsll', 'popcountll']
for builtin in sabnatkBuiltins:
  if not conf.CheckDeclaration('__builtin_%s' % builtin):
    print('ERROR: __builtin_%s is required by SABNAtk')
    Exit(1)
# Check for aligned malloc declaration
if conf.CheckCHeader('mm_malloc.h'):
  if not conf.CheckDeclaration('_mm_malloc', '#include <mm_malloc.h>'):
    print('ERROR: _mm_malloc is required by SABNAtk')
    Exit(1)
else:
  print('ERROR: mm_malloc.h not found')
  Exit(1)
conf.env.Append(CXXFLAGS='-msse4')
if not conf.CheckCXX():
  conf.env.Replace(CXXFLAGS=cppFlags)
  print('WARNING: -msse4 not supported; SABNAtk performance will be affected')
# Check for boost header location
if not conf.CheckCXXHeader('boost/version.hpp'):
  Exit(1)
# Check for all the libraries
for lib in allLibs:
  if not conf.CheckLib(lib, language='C++'):
    Exit(1)
# Use the modified environment
env = conf.Finish()

env.targetName = targetName
env.testName = testName
env.topDir = topDir
env.testDir = os.path.join(topDir, 'test')

buildDir = os.path.join('builds', buildDir)

if not SConscript('SConscript', exports='env', src_dir='.', variant_dir=buildDir, duplicate=0):
  print('Executable was not built')
  Exit(1)
if ARGUMENTS.get('TEST', 1) not in ('0', 0):
  if not SConscript(os.path.join('test', 'SConscript'), exports='env', src_dir='test', variant_dir=os.path.join(buildDir, 'test'), duplicate=0):
    print('Test suite was not built')
