import stat
import os

from typing import List, Union, Optional
from os.path import join as opj
from itertools import chain
from collections.abc import Iterable

from .crs import CRSResults
from ...core.settings import Settings
from ...core.basejob import SingleJob
from ...core.errors import JobError, ResultsError, FileError
from ...mol.molecule import Molecule

try:
    from rdkit.Chem import MolToSmiles, RemoveHs
    from ..molecule.rdkit import from_smiles
    _RD_EX = False
except ImportError as ex:
    _RD_EX = ex

__all__ = ['UnifacJob', 'UnifacResults']


class UnifacResults(CRSResults):
    """The :class:`Results` subclass assigned to :class:`UnifacJob`."""

    def recreate_molecule(self) -> List[Molecule]:
        """Reconstruct and return list with all input molecules.

        Molecules are extracted from the SMILES string(s) in the .run file.
        """
        if _RD_EX:
            err = f"SMILES to Molecule conversion requires the 'rdkit' package: {_RD_EX}"
            raise ResultsError(err)

        s = self.recreate_settings()
        if isinstance(s.input.smiles, list):
            return [from_smiles(smiles) for smiles in s.input.smiles]
        return [from_smiles(s.input.smiles)]

    def recreate_settings(self) -> Settings:
        """Reconstruct a :class:`Settings` instance from the .run file."""
        def _read_runfile(runfile: str) -> Optional[List[str]]:
            """Extract the command line input from the runfile"""
            ret = None
            with open(runfile, 'r') as f:
                for i in f:
                    if '"$ADFBIN"/unifac' in i:
                        ret = i.split()
                        arg_list[-1] = arg_list[-1].rstrip('\\')
                    else:
                        continue

                    while i.endswith('\\'):  # THe input might be spread over multiple lines
                        i = next(f)
                        ret += i.split().rstrip('\\')
                    del ret[0]  # Delete ``"$ADFBIN"/unifac``
                    break
            return ret

        def _runfile2settings(arg_list: List[str]) -> Settings:
            """Parse the content of the extracted .run file."""
            iterator = iter(arg_list)
            key = None
            s = Settings()
            for value in iterator:
                # Identify potential keys
                if value.startswith('-'):
                    key = value.strip('-').lower()
                    continue

                #  If possible, convert value into a float or integer
                value = self._str_to_number(value)

                # Create a new key/value pair or update an old key/value pair with a list of values
                if key in s.input:
                    try:
                        s.input[key].append(value)
                    except AttributeError:
                        s.input[key] = [s.input.pop(key), value]
                else:
                    s.input[key] = value
            return s

        # Read the .run file
        runfile = self['$JN.run']
        arg_list = _read_runfile(runfile)
        if arg_list is None:
            raise FileError(f"recreate_settings: Failed to parse the content of '{runfile}'")

        # Reconstruct and return the job settings
        return _runfile2settings(arg_list)

    @staticmethod
    def _str_to_number(value: str) -> Union[str, int, float]:
        """Attempt to convert *value* into a :class:`float` or :class:`int`."""
        if '.' in value:
            try:
                return float(value)
            except ValueError:
                return value
        else:
            try:
                return int(value)
            except ValueError:
                return value


class UnifacJob(SingleJob):
    """A :class:`Job` subclass for interfacing with ADFs' implementation of UNIFAC_.

    .. _UNIFAC: https://www.scm.com/doc/COSMO-RS/The_UNIFAC_program.html

    """

    _result_type = UnifacResults

    def __init__(self, **kwargs: object) -> None:
        super().__init__(**kwargs)

        # If supplied, convert self.molecule into a SMILES string
        if self.molecule and _RD_EX:
            err = f"Molecule to SMILES conversion requires the 'rdkit' package: {_RD_EX}"
            raise JobError(err)

        elif self.molecule:
            mol_list = [self.molecule] if isinstance(self.molecule, Molecule) else self.molecule
            k1 = self.settings.input.find_case('smiles')
            k2 = self.settings.input.find_case('-smiles')
            smiles = k1 if k2 not in self.settings.input else k2
            self.settings.input[smiles] = [MolToSmiles(RemoveHs(mol)) for mol in mol_list]

    def _get_ready(self) -> None:
        """Create the runfile."""
        runfile = opj(self.path, self._filename('run'))
        with open(runfile, 'w') as run:
            run.write(self.full_runscript())
        os.chmod(runfile, os.stat(runfile).st_mode | stat.S_IEXEC)

    def get_input(self) -> None: return None
    def hash_input(self) -> str: return self.hash_runscript()

    def get_runscript(self) -> str:
        """Run a MACTH runscript."""
        iterator = self.settings.input.items()
        kwargs = {self._sanitize_key(k): self._sanitize_value(v) for k, v in iterator}
        kwargs_iter = chain.from_iterable(kwargs.items())
        args = ' '.join(i for i in kwargs_iter)
        return f'"$ADFBIN"/unifac {args}'

    """###################################### New methods ######################################"""

    @staticmethod
    def _sanitize_key(key: str) -> str:
        """Lower *key* and prepended it with the ``"-"`` character."""
        try:
            ret = key if key.startswith('-') else f'-{key}'
        except AttributeError:  # *key* is not a string
            raise JobError(f"Key of invalid type encountered in Job settings: {repr(key)}")
        return ret.lower()

    @staticmethod
    def _sanitize_value(value: object) -> str:
        """Convert *value* into a string; join its elements if it's an iterable."""
        if isinstance(value, str) or not isinstance(value, Iterable):
            return repr(value)
        return ' '.join(repr(i) for i in value)
