##
# @file SConscript
# @brief Provides functionality for building the executable.
# @author Ankit Srivastava <asrivast@gatech.edu>
#
# Copyright 2020 Georgia Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
Import('env')

boostLibs = [
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
  if lib in buildEnv.get('LIBS', []):
    continue
  if not conf.CheckLib(lib, language='C++'):
    Return('built')
buildEnv = conf.Finish()

buildEnv.Program(target=buildEnv.targetName, source=srcFiles)

buildEnv.Install(buildEnv.topDir, buildEnv.targetName)
built = True
Return('built')
