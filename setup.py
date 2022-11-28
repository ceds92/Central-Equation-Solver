#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = '0.1.0'

setup(
    name='Central-Equation-Solver',
    version=version,
    author='Julian Ceddia',
    author_email='jdceddia@gmail.com',
    description='Solve the 2D TISE in a periodic potential',
    long_description=long_description,
    url='https://github.com/ceds92/Central-Equation-Solver',
    project_urls = {
        "Bug Tracker": "https://github.com/ceds92/Central-Equation-Solver/issues"
    },
    license='MIT',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'customtkinter','matplotlib_scalebar','scipy'],
)
