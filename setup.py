import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for CPAT  -- Coding Potential Assessment Tool
"""

def main():
	setup(  name = "CPAT",
			version = "3.0.0",
			py_modules = [ 'psyco_full' ],
			python_requires='>=3.5',
			packages = find_packages( 'lib' ),
			package_dir = { '': 'lib' },
			package_data = { '': ['*.ps'] },
			scripts = glob.glob( "bin/*.py"),
			ext_modules = [],
			test_suite = 'nose.collector',
			setup_requires = ['nose>=0.10.4'],
			author = "Liguo Wang, Jung Hyun Park",
			author_email ="wangliguo78@gmail.com, hj_park@pitt.edu",
			platforms = ['Linux','MacOS'],
			requires = ['cython (>=0.17)'],
			install_requires = ['numpy','pysam'], 
			description = "CPAT (Coding Potential Assessment Tool)",
			long_description = "CPAT is an alignment-free method to predict RNA coding potential using four sequence features.",
			license='GNU General Public License',
			url = "http://rna-cpat.sourceforge.net/",
			zip_safe = False,
			dependency_links = [],
			classifiers=[
				'Development Status :: 5 - Production/Stable',
				'Environment :: Console',
				'Intended Audience :: Science/Research',
				'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: POSIX',
				'Programming Language :: Python',
				'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
			
			keywords='RNA coding potential prediction, lncRNA, lincRNA, logsitic regression',
             )


if __name__ == "__main__":
	main()
