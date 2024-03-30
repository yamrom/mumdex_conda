#! /usr/bin/env python

from setuptools import setup, Extension, find_packages

# Define the extension module
extension = Extension('mumdex._mumdex',
                       sources=['src/python_mumdex.cpp',
                                'src/files.cpp',
                                'src/utility.cpp'],
                       include_dirs=['src'],
                       extra_compile_args=['-std=c++11',
                                           '-pthread',
                                           '-Isrc'])

setup(
    name='mumdex',
    version='0.1.0',
    description="MUMdex python module",
    license="MIT",
    author="Boris Yamrom",
    author_email='yamrom.boris@gmail.com',
    url='https://github.com/iossifovlab/mumdex',
    ext_modules=[extension],
    scripts=['scripts/test_python_mumdex.sh',
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
             'scripts/bridge_tester.py'],
    packages=find_packages(),
    requires=['cmake', 'numpy', 'sys', 'os']
)
