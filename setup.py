from distutils.core import setup
import dimepy
import sys

setup(
    name='dimepy',
    version='0.0.7',
    packages=['dimepy'],
    url='http://www.github.com/KeironO/dimepy',
    license='GPLv2',
    platforms=['Windows', 'UNIX'],
    install_requires=open('requirements.txt').read().splitlines(),
    long_description=open('README.md').read(),
    author='Keiron O\'Shea',
    author_email = 'keo7@aber.ac.uk',
    description = 'Python package for the high-thoroughput nontargeted metabolite fingerprinting of nominal mass direct injection mass spectrometry.',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers'
    ]
)
