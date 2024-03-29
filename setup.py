#! /usr/bin/env python

from setuptools import setup, Extension



# Define the extension module
extension = Extension('mumdex._mumdex',
                       sources=['src/python_mumdex.cpp', 'src/files.cpp','src/utility.cpp'],
                       include_dirs=['src'],
                       extra_compile_args=['-std=c++11', '-pthread', '-Isrc'])
#,
#                       extra_link_flags=['-lbin/libfiles.a','-lbin/utility.a'] #)

setup(
    name='mumdex',
    version='0.1.0',
    description='Your package description',
    ext_modules=[extension],
    data_files=[('bin',
                   ['scripts/test_python_mumdex.sh',
                    'scripts/bridge_finder.py',
                    'scripts/bridge_info.py',
                    'scripts/candidate_finder.py',
                    'scripts/chromosome_bridges.py',
                    'scripts/load_counts.py',
                    'scripts/mumdex2txt.py',
                    'scripts/show_mums.py',
                    'scripts/test_python_mumdex.py',
                    'scripts/mapper.py',
                    'scripts/sample_bridge_finder.sh',
                    'scripts/bridge_tester.py'])],
      requires=['numpy', 'sys', 'os']

)



