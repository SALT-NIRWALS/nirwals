from setuptools import setup, find_namespace_packages

# with open("README.rst", "rb") as f: #
#     long_descr = f.read().decode("utf-8")

setup(
    name = "nirwals",
    package_dir = {'': 'src'},
    packages = ['nirwals'],
    version = '{{VERSION_PLACEHOLDER}}',
    description = "A pipeline for instrument signature removal for the SALT/NIRWALS instrument",
    # # long_description = long_descr,
    author = "Ralf Kotulla",
    author_email = "kotulla@wisc.edu",
    url = "https://github.com/SALT-NIRWALS/nirwals",

    scripts=['src/nirwals_reduce.py',
             'src/nirwals_makedark.py',
             'src/nirwals_makemasterdark.py',
             'src/nirwals_fit_nonlinearity.py',
             ],

    install_requires=[
        'astropy',
        'matplotlib',
        'multiparlog',
        'numpy',
        'pandas',
        'scipy',
        'pyvo',
        'watchdog',
        ],

    )