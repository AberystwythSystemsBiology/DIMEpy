try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

try:
    import pypandoc
    long_description = pypandoc.convert("README.md", "rst")
except (IOError, ImportError):
    long_description=open("README.md").read()



setup(
    name='dimepy',
    version='0.0.9',
    url='http://www.github.com/KeironO/dimepy',
    license='GPLv2',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=open("requirements.txt").read().splitlines(),
    platforms=['Windows', 'UNIX'],
    long_description=long_description,
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
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers'
    ]
)
