"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
import os

here = os.path.abspath(os.path.dirname(__file__))
bindir = os.path.join(here, "bin/")

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
	long_description = f.read()

setup(name='pygpseq',
	version='3.3.4',
	description='A GPSeq image analysis package',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ggirelli/gpseq-img-py',
	author='Gabriele Girelli',
	author_email='gabriele.girelli@scilifelab.se',
	license='MIT',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3 :: Only',
	],
	keywords='microscopy image analysis bioimaging biology cell DNA',
	packages=find_packages(),
	install_requires=[
		'ggc>=0.0.3',
		'czifile>=0.1.5',
		'cython',
		'jinja2==2.10',
		'joblib==0.11',
		'matplotlib==2.2.2',
		'nd2reader>=3.1.0',
		'numpy>=1.14.2',
		'pandas>=0.22.0',
		'scipy>=1.0.0',
		'scikit-image==0.14.0',
		'seaborn>=0.9.0',
		'tifffile>=0.15.1',
		'tqdm>=4.23.4',
		'weasyprint>=0.42.2'],
	scripts=[os.path.join(bindir, fp) for fp in os.listdir(bindir)],
	test_suite="nose.collector",
	tests_require=["nose"],
)
