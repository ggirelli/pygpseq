"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

setup(name='pygpseq',
	version='1.0.0',
	description='A GPSeq image analysis package',
	long_description=long_description,
	url='https://github.com/ggirelli/gpseq-img-py',
	author='Gabriele Girelli',
	author_email='gabriele.girelli@scilifelab.se',
	license='MIT',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3'
	],
	keywords='microscopy image analysis bioimaging biology cell DNA',
	packages=["pygpseq"],
	install_requires=['jinja2', 'joblib', 'matplotlib', 'numpy', 'pandas',
	'scipy', 'scikit-image', 'tifffile', 'weasyprint'],
	scripts=["bin/gpseq_anim"],
	test_suite="nose.collector",
	tests_require=["nose"],
)
