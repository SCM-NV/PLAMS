from os.path import join as opj

from setuptools import find_packages, setup

with open(opj("version.py")) as f:
    exec(f.read())

packages = ["scm.plams"] + ["scm.plams." + i for i in find_packages(".")]

description = "PLAMS is a library providing powerful, flexible and easily extendable Python interface to molecular modeling programs. It takes care of input preparation, job execution, file management and output data extraction. It helps with building advanced data workflows that can be executed in parallel, either locally or by submitting to a resource manager queue.\n\nPlease check the project's GitHub page for more information: https://github.com/SCM-NV/PLAMS \n\nPLAMS is an Open Source project supported by `Software for Chemistry & Materials B.V. <https://www.scm.com>`_"

setup(
    name="plams",
    version=__version__,
    author="Software for Chemistry and Materials",
    author_email="info@scm.com",
    url="https://www.scm.com/doc/plams/",
    download_url="https://github.com/SCM-NV/PLAMS/zipball/release",
    license="LGPLv3",
    description="Python Library for Automating Molecular Simulations",
    long_description=description,
    classifiers=[
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords=["molecular modeling", "computational chemistry", "workflow", "python interface"],
    python_requires=">=3.6",
    install_requires=["dill>=0.2.4", "numpy<2", "scipy", "natsort"],
    extras_require={"chem_tools": ["ase<3.23", "rdkit"], "results": ["networkx<3"]},
    packages=packages,
    package_dir={"scm.plams": "."},
    package_data={"scm.plams": [".flake8", "plams_defaults", "examples/*", "unit_tests/*", "unit_tests/xyz/*"]},
    scripts=[opj("scripts", "plams")],
)
