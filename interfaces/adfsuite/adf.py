import warnings

import numpy as np

from scm.plams.core.errors import ResultsError
from scm.plams.core.functions import config, log
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.scmjob import SCMJob, SCMResults
from scm.plams.tools.periodic_table import PT
from scm.plams.tools.units import Units

__all__ = ["ADFJob", "ADFResults"]


class ADFResults(SCMResults):
    """A specialized |SCMResults| subclass for accessing the results of |ADFJob|."""

    _kfext = ".t21"
    _rename_map = {"TAPE{}".format(i): "$JN.t{}".format(i) for i in range(10, 100)}

    def get_properties(self):
        """get_properties()
        Return a dictionary with all the entries from ``Properties`` section in the main KF file (``$JN.t21``).
        """
        ret = {}
        for sec, var in self._kf:
            if sec == "Properties":
                ret[var] = self.readkf(sec, var)
        return ret

    def get_main_molecule(self):
        """get_main_molecule()
        Return a |Molecule| instance based on the ``Geometry`` section in the main KF file (``$JN.t21``).

        For runs with multiple geometries (geometry optimization, transition state search, intrinsic reaction coordinate) this is the **final** geometry. In such a case, to access the initial (or any intermediate) coordinates please use :meth:`get_input_molecule` or extract coordinates from section ``History``, variables ``xyz 1``, ``xyz 2`` and so on. Mind the fact that all coordinates written by ADF to ``History`` section are in bohr and internal atom order::

            mol = results.get_molecule(section='History', variable='xyz 1', unit='bohr', internal=True)
        """
        return self.get_molecule(section="Geometry", variable="xyz InputOrder", unit="bohr")

    def get_input_molecule(self):
        """get_input_molecule()
        Return a |Molecule| instance with initial coordinates.

        All data used by this method is taken from ``$JN.t21`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        if ("History", "nr of geometries") in self._kf:
            return self.get_molecule(section="History", variable="xyz 1", unit="bohr", internal=True)
        return self.get_main_molecule()

    def get_energy(self, unit="au"):
        """get_energy(unit='au')
        Return final bond energy, expressed in *unit*.
        """
        return self._get_single_value("Energy", "Bond Energy", unit)

    def get_dipole_vector(self, unit="au"):
        """get_dipole_vector(unit='au')
        Return the dipole vector, expressed in *unit*.
        """
        prop = self.get_properties()
        if "Dipole" in prop:
            return Units.convert(prop["Dipole"], "au", unit)
        raise ResultsError("'Dipole' not present in 'Properties' section of {}".format(self._kfpath()))

    def get_gradients(self, energy_unit="au", dist_unit="bohr", eUnit=None, lUnit=None):
        """get_gradients(energy_unit='au', dist_unit='bohr')
        Return the cartesian gradients from the 'Gradients_InputOrder' field of the 'GeoOpt' Section in the kf-file, expressed in given units. Returned value is a numpy array with shape (nAtoms,3).
        """
        if eUnit:
            log(
                "Deprecated Keyword eUnit used in ADFResults.get_gradients, update your script! Overwriting energy_unit with the given argument.",
                1,
            )
            warnings.warn("eUnit is deprecated, use energy_unit instead.", category=DeprecationWarning)
            energy_unit = eUnit
        if lUnit:
            log(
                "Deprecated Keyword lUnit used in ADFResults.get_gradients, update your script! Overwriting energy_unit with the given argument.",
                1,
            )
            warnings.warn("lUnit is deprecated, use dist_unit instead.", category=DeprecationWarning)
            dist_unit = lUnit
        gradients = np.array(self.readkf("GeoOpt", "Gradients_InputOrder"))
        gradients.shape = (-1, 3)
        gradients *= Units.conversion_ratio("au", energy_unit) / Units.conversion_ratio("bohr", dist_unit)
        return gradients

    def _extract_hessian(self, section, variable, internal_order):
        """_extract_hessian(section, variable, internal_order)
        Extract Hessian from *section*/*variable* of the TAPE21 file. Reorder from internal to input order, if *internal_order* is ``True``.
        """
        hess_int = np.array(self.readkf(section, variable))
        n = int((len(hess_int) / 9 + 1) ** 0.5)
        hess_int.shape = (3 * n, 3 * n)
        if internal_order:
            hess_inp = np.zeros(hess_int.shape)
            mapping = self._int2inp()
            for i in range(n):
                for j in range(n):
                    ii, jj = mapping[i] - 1, mapping[j] - 1
                    hess_inp[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = hess_int[3 * ii : 3 * ii + 3, 3 * jj : 3 * jj + 3]
            return hess_inp
        return hess_int

    def get_hessian(self):
        """get_hessian()
        Try extracting Hessian, either analytical or numerical, whichever is present in the TAPE21 file, in the input order. Returned value is a square numpy array of size 3*nAtoms.
        """
        if ("Hessian", "Analytical Hessian") in self._kf:
            return self._extract_hessian("Hessian", "Analytical Hessian", True)
        if ("Freq", "Hessian_complete") in self._kf:
            return self._extract_hessian("Freq", "Hessian_complete", True)
        raise ResultsError("auto_hessian: Hessian does not seem to be present in t21 file.")

    def get_energy_decomposition(self, unit="au"):
        """get_energy(unit='au')
        Return a dictionary with energy decomposition terms, expressed in *unit*.

        The following keys are present in the returned dictionary: ``Electrostatic``, ``Kinetic``, ``Coulomb``, ``XC``. The sum of all the values is equal to the value returned by :meth:`get_energy`.
        Note that additional contributions might be included, those are up to now: ``Dispersion``.
        """
        ret = {}
        ret["Electrostatic"] = self._get_single_value("Energy", "Electrostatic Energy", unit)
        ret["Kinetic"] = self._get_single_value("Energy", "Kinetic Energy", unit)
        ret["Coulomb"] = self._get_single_value("Energy", "Elstat Interaction", unit)
        ret["XC"] = self._get_single_value("Energy", "XC Energy", unit)
        if ("Energy", "Dispersion Energy") in self._kf:
            ret["Dispersion"] = self._get_single_value("Energy", "Dispersion Energy", unit)
        return ret

    def get_frequencies(self, unit="cm^-1"):
        """get_frequencies(unit='cm^-1')
        Return a numpy array of vibrational frequencies, expressed in *unit*.
        """
        freqs = np.array(self.readkf("Freq", "Frequencies"))
        return freqs * Units.conversion_ratio("cm^-1", unit)

    def get_timings(self):
        """get_timings()

        Return a dictionary with timing statistics of the job execution. Returned dictionary contains keys ``cpu``, ``system`` and ``elapsed``. The values are corresponding timings, expressed in seconds.
        """
        last = self.grep_output(" Total Used : ")[-1].split()
        ret = {}
        ret["elapsed"] = float(last[-1])
        ret["system"] = float(last[-3])
        ret["cpu"] = float(last[-5])
        return ret

    def _atomic_numbers_input_order(self):
        """_atomic_numbers_input_order()
        Return a list of atomic numbers, in the input order.
        """
        n = self.readkf("Geometry", "nr of atoms")
        tmp = self.readkf("Geometry", "atomtype").split()
        atomtypes = {i + 1: PT.get_atomic_number(tmp[i]) for i in range(len(tmp))}
        atomtype_idx = self.readkf("Geometry", "fragment and atomtype index")[-n:]
        atnums = [atomtypes[i] for i in atomtype_idx]
        return self.to_input_order(atnums)

    def _int2inp(self):
        """_int2inp()
        Get mapping from the internal atom order to the input atom order.
        """
        aoi = self.readkf("Geometry", "atom order index")
        n = len(aoi) // 2
        return aoi[:n]

    def recreate_molecule(self):
        """Recreate the input molecule for the corresponding job based on files present in the job folder. This method is used by |load_external|."""
        if self._kfpresent():
            return self.get_input_molecule()
        return None

    def recreate_settings(self):
        """Recreate the input |Settings| instance for the corresponding job based on files present in the job folder. This method is used by |load_external|."""
        if self._kfpresent():
            if ("General", "Input") in self._kf:
                tmp = self.readkf("General", "Input")
                user_input = "\n".join([tmp[i : 160 + i].rstrip() for i in range(0, len(tmp), 160)])
            else:
                user_input = self.readkf("General", "user input")
            try:
                from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade

                inp = InputParserFacade().to_settings("adf", user_input)
            except:
                log("Failed to recreate input settings from {}".format(self._kf.path), 5)
                return None

            s = Settings()
            s.input = inp
            del s.input.atoms
            s.soft_update(config.job)
            return s
        return None


class ADFJob(SCMJob):
    _result_type = ADFResults
    _command = "adf"

    def _serialize_mol(self):
        for i, atom in enumerate(self.molecule):
            smb = self._atom_symbol(atom)
            suffix = ""
            if "adf" in atom.properties and "fragment" in atom.properties.adf:
                suffix += "f={fragment} "
            if "adf" in atom.properties and "block" in atom.properties.adf:
                suffix += "b={block}"

            self.settings.input.atoms["_" + str(i + 1)] = ("{:>5}".format(i + 1)) + atom.str(
                symbol=smb, suffix=suffix, suffix_dict=atom.properties.adf
            )

    def _remove_mol(self):
        if "atoms" in self.settings.input:
            del self.settings.input.atoms
