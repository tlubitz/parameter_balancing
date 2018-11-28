"""
Setup file for parameter balancing standalone version.
The standalone version can be integrated into existing workflows of Systems Biology.
See: https://bitbucket.org/tlubitz/pb

"""

from setuptools import setup,find_packages

setup(name='pbalancing',version='2.0.37',description='Parameter Balancing for Kinetic Metabolic Models',long_description='Parameter Balancing - A Systems Biology tool to determine thermodynamically consistent parameter sets for kinetic models of cell metabolism',url='http://www.parameterbalancing.net',author='Timo Lubitz',author_email='timo.lubitz@gmail.com',license='MIT',classifiers=['Development Status :: 4 - Beta','Intended Audience :: Developers', 'Topic :: Software Development', 'Programming Language :: Python :: 3.0'],keywords='modelling systems biology kinetic parameter estimation thermodynamics',packages=find_packages(),install_requires=['python-libsbml','tablib'],python_requires='>=3.0',include_package_data=True)
