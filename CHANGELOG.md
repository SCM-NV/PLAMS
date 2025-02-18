# Changelog
Notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

The versioning for this project aligns with that of AMS. 
The format for this is the original year of the main release as a prefix, and an incremental postfix for each sub-release, starting at 101.
In general, subsequent sub-releases after the main release contain only bug fixes.
For example, 2024.101 is the major release of 2024, and 2024.102 is the first bugfix release.

This changelog is effective from the 2025 releases.

## [Unreleased]

### Added
* Methods `get_system`, `get_input_system` and `get_main_system` to `AMSResults`, which return an AMS `ChemicalSystem` instead of a PLAMS `Molecule` 
* `AMSJob` can accept an AMS `ChemicalSystem` instead of a PLAMS `Molecule` as an input system
* Specific `ConfigSettings` and related settings classes with explicitly defined fields
* Support for work functions: `AMSResults.get_work_function_results` and `plot_work_function`
* New `packmol_around` function for packing in non-orthorhombic boxes.
* New `plot_grid_molecules` function for plotting with rdkit multiple molecules.
* `Molecule.delete_atoms` method to delete multiple atoms with partial success 
* Examples on `MoleculeFormats` and `MoleculeTools`
* Examples on `Logging`
* Script `generate_example.sh` to generate documentation pages from notebook examples
* GitHub workflows for CI and publishing to PyPI
* Build using `pyproject.toml`, addition of extras groups to install optional dependencies
* Logging of job summaries to CSV logfile
* Logging of AMS job error messages to stdout and logfile on job failure
* Method `get_errormsg` enforced on the `Job` base class, with a default implementation
* Added an interface to the Serenity program through methods such as `SerenityJob`, `SerenityResults` and `SerenitySettings`
* Methods `Settings.nested_keys`, `Settings.block_keys` for accessing nested keys and `Settings.contains_nested` and `Settings.pop_nested` for checking if nested keys exist in a settings object and popping them
* Method `Settings.compare` added to compare and contrast items in two settings objects
* Method `readcoskf` added to `Molecule` class, enabling the reading of COSKF file

### Changed
* Functions for optional packages (e.g. RDKit, ASE) are available even when these packages are not installed, but will raise an `MissingOptionalPackageError` when called
* `AMSResults.get_main_ase_atoms` also includes atomic charges
* Global `config` is initialized with a `ConfigSettings` instead of loading from the standard `plams_defaults` file
* `init` and `finish` functions are now optional
* `Job.status` is a `JobStatus` string enum
* Supercell and RDKit properties are no longer serialized to AMS input
* Restructuring of examples and conversion of various examples to notebooks
* Support for `networkx>=3` and `ase>=3.23`
* Use standard library logger for `log` function
* Make `Job` class inherit from `ABC` and mark abstract methods 
* Exceptions raised in `prerun` and `postrun` will always be caught and populate error message
* `Settings.get_nested` takes a default argument which is returned if the nested key is not present in the settings instance

### Fixed
* `Molecule.properties.charge` is a numeric instead of string type when loading molecule from a file
* `Molecule.delete_all_bonds` removes the reference molecule from the removed bond instances
* `SingleJob.load` returns the correctly loaded job
* `AMSJob.check` handles a `NoneType` status, returning `False`
* `MultiJob.run` locking resolved when errors raised within `prerun` and `postrun` methods
* `Molecule.add_hatoms` to use bonding information if available when adding new hydrogen atoms

### Deprecated
* `plams` launch script is deprecated in favor of simply running with `amspython`

### Removed
* Legacy `BANDJob`, `DFTBJob`, `UFFJob`, `MOPACJob`, `ReaxFFJob`, `CSHessianADFJob` and `ADFJob` have been removed
* Exception classes `AMSPipeDecodeError`, `AMSPipeError`, `AMSPipeInvalidArgumentError`, `AMSPipeLogicError`, `AMSPipeRuntimeError`, `AMSPipeUnknownArgumentError`, `AMSPipeUnknownMethodError`, `AMSPipeUnknownVersionError`, were moved from scm.plams to scm.amspipe.




