import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "fvm",
    version = "0.1",
    author = "Daniel J Farrell",
    author_email = "daniel.james.farrell@gmail.com",
    description = ("Examples of solving conservation problems using the finite volume method"),
    license = "BSD",
    keywords = "finite volume numerical partial differential equations",
    url = "http://github.com/danieljfarrell/FVM",
    packages=['fvm'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: BSD License",
    ],
)