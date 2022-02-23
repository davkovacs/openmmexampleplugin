from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
ACEplugin_header_dir = '@ACEPLUGIN_HEADER_DIR@'
ACEplugin_library_dir = '@ACEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_ACEplugin',
                      sources=['ACEPluginWrapper.cpp'],
                      libraries=['OpenMM', 'ACEPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), ACEplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), ACEplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='ACEplugin',
      version='1.0',
      py_modules=['ACEplugin'],
      ext_modules=[extension],
     )
