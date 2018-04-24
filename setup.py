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
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

setup(name='pygpseq',
	version='2.1.3',
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
	install_requires=['jinja2', 'joblib', 'matplotlib', 'nd2reader', 'numpy',
	'pandas', 'scipy', 'scikit-image', 'tifffile', 'tqdm', 'weasyprint'],
	scripts=[os.path.join(bindir, fp) for fp in os.listdir(bindir)],
	test_suite="nose.collector",
	tests_require=["nose"],
)
