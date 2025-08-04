import shutil
from os import symlink
from os.path import basename
from os.path import join as opj
from os.path import relpath

import numpy as np

from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import JobError
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units

__all__ = ["ORCAJob", "ORCAResults"]


class ORCAResults(Results):
    """A class for ORCA results."""

    def get_runtime(self):
        """Return runtime in seconds from output."""
        runtimeString = self.grep_output("TOTAL RUN TIME")[-1]
        runtimeList = [int(x) for x in runtimeString.split()[3::2]]
        assert len(runtimeList) == 5  # days hours minutes seconds milliseconds
        runtime = sum([t * m for t, m in zip(runtimeList, [24 * 60 * 60, 60 * 60, 60, 1, 0.001])])
        return runtime

    def get_timings(self):
        """Return timings section as dictionary. Units are seconds."""
        timingsSection = self.get_output_chunk("Timings for individual modules:", end="****")
        ret = {}
        for line in timingsSection:
            line = line.strip()
            if not line:
                continue
            name, data = line.split("...")
            ret[name.strip()] = float(data.split()[0])
        return ret

    def check(self):
        """Returns true if ORCA TERMINATED NORMALLY is in the output"""
        return bool(self.grep_output("ORCA TERMINATED NORMALLY"))

    def check_go_conv(self):
        """Returns true if THE OPTIMIZATION HAS CONVERGED is in the output"""
        return bool(self.grep_output("THE OPTIMIZATION HAS CONVERGED"))

    def check_scf_conv(self):
        """Returns true if SCF NOT CONVERGED AFTER is NOT in the output"""
        return not bool(self.grep_output("SCF NOT CONVERGED AFTER"))

    def _get_energy_correct_unit(self, string, unit="Eh"):
        # if there are multiple numbers in the string, get the one
        # labeled with unit otherwise get the last number
        stringList = string.split()
        if unit in string:
            index = stringList.index(unit)
            return float(stringList[index - 1])
        else:
            return float(stringList[-1])

    def _get_energy_type(self, search="FINAL SINGLE POINT ENERGY", index=-1, unit="a.u.", filterDotDotDot=True):
        # hacky way of getting rid of some entries: first get all, remove
        # those that are '... done'
        s = self.grep_output(search)
        if filterDotDotDot:
            s = [item for item in s if "..." not in item]
        s = s[index]
        if not isinstance(index, slice):
            return Units.convert(self._get_energy_correct_unit(s), "a.u.", unit)
        else:
            return [Units.convert(self._get_energy_correct_unit(x), "a.u.", unit) for x in s]

    def get_scf_iterations(self, index=-1):
        """Returns Number of SCF Iterations from the Output File.

        Set ``index`` to choose the n-th occurence, *e.g.* to choose
        an certain step. Also supports slices.
        Defaults to the last occurence.
        """
        s = self.grep_output("SCF CONVERGED AFTER")
        n = [int(x.split()[-3]) for x in s]
        return n[index]

    def get_energy(self, index=-1, unit="a.u."):
        """Returns 'FINAL SINGLE POINT ENERGY:' from the output file.

        Set ``index`` to choose the n-th occurence of the total energy
        in the output, *e.g.* to choose an certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("FINAL SINGLE POINT ENERGY", index=index, unit=unit)

    def get_dispersion(self, index=-1, unit="a.u."):
        """Returns 'Dispersion correction' from the output file.

        Set ``index`` to choose the n-th occurence of the dispersion energy
        in the output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Dispersion correction", index=index, unit=unit)

    def get_zpe(self, index=-1, unit="a.u."):
        """Returns 'Non-thermal (ZPE) correction' from the output file.

        Set ``index`` to choose the n-th occurence of the ZPE in the output,
        *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Non-thermal (ZPE) correction", index=index, unit=unit)

    def get_inner_energy(self, index=-1, unit="a.u."):
        """Returns 'Total thermal energy' from the output file.

        Set ``index`` to choose the n-th occurence of the inner energy in
        the output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Total thermal energy", index=index, unit=unit)

    def get_enthalpy(self, index=-1, unit="a.u."):
        """Returns 'Total Enthalpy' from the output file.

        Set ``index`` to choose the n-th occurence of the enthalpy in the
        output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Total Enthalpy", index=index, unit=unit, filterDotDotDot=False)

    def get_entropy(self, index=-1, unit="a.u."):
        """Returns 'Final entropy term' from the output file.

        Set ``index`` to choose the n-th occurence of the entropy in the
        output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Final entropy term", index=index, unit=unit, filterDotDotDot=False)

    def get_gibbs_free_energy(self, index=-1, unit="a.u."):
        """Returns 'Final Gibbs free energy' from the output file.

        Set ``index`` to choose the n-th occurence of the Gibbs free energy in
        the output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        return self._get_energy_type("Final Gibbs free energy", index=index, unit=unit, filterDotDotDot=False)

    def get_electrons(self, index=-1, spin_resolved=False):
        """Get Electron count from Output.

        Set ``spins`` to ``True`` to get a tuple of alpha and beta electrons
        instead of totals.  Set ``index`` to choose the n-th occurcence in the
        output, *e.g.* to choose a certain step.
        Also supports slices.
        Defaults to the last occurence.
        """
        if spin_resolved:
            alpha = [float(s.split()[-2]) for s in self.grep_output("N(Alpha)")][index]
            beta = [float(s.split()[-2]) for s in self.grep_output("N(Beta)")][index]
            if isinstance(alpha, float):
                alpha = [alpha]
                beta = [beta]
            ret = [(a, b) for a, b in zip(alpha, beta)]
        else:
            ret = [float(s.split()[-2]) for s in self.grep_output("N(Total)")][index]
        return ret

    def get_gradients(self, match=0, energy_unit="a.u.", dist_unit="bohr"):
        """Returns list of ndarrays with forces from the output
        (there the unit is a.u./bohr).

        ``match`` is passed to :meth:`~Results.get_output_chunk`,
        defaults to 0.
        """
        conv = Units.conversion_ratio("a.u.", energy_unit) / Units.conversion_ratio("bohr", dist_unit)
        searchBegin = "CARTESIAN GRADIENT"
        searchEnd = "Difference to translation invariance:"
        block = self.get_output_chunk(begin=searchBegin, end=searchEnd, match=match)
        ret = []
        for line in block:
            line = line.strip().split()
            if len(line) == 1:  # ---- lines means new set of forces
                ret.append([])
                continue
            if not line:  # ignore empty line
                continue
            ret[-1].append(line[-3:])
        ret = [np.array(item, dtype=float) * conv for item in ret]
        if len(ret) == 1:
            return ret[0]
        else:
            return ret

    def get_dipole_vector(self, index=-1, unit="a.u."):
        """Get the Dipole Vector
        Returns the dipole vector, expressed in *unit*.
        """
        lines = self.grep_output("Total Dipole Moment")
        conv = Units.conversion_ratio("a.u.", unit)
        vectors = [np.array(line.split()[-3:], dtype=float) * conv for line in lines]
        return vectors[index]

    def get_dipole(self, **kwargs):
        """Get Hirshfeld Analysis from Output

        Uses :meth:`~ORCAResults.get_dipole_vector` to calculate the
        total dipole.
        All options are passed on.
        """
        vec = self.get_dipole_vector(**kwargs)
        if isinstance(vec, np.ndarray):
            vec = [vec]
        vec = np.linalg.norm(vec, axis=1)
        if len(vec) == 1:
            return vec[0]
        else:
            return vec

    def get_atomic_charges(self, method="mulliken", match=0):
        """Get Atomic Charges from Output

        :parameter method: Can be any one that is available in the output, e.g. mulliken or loewdin.
        :parameter match: Select occurence in the output to use. E.g. when running multiple structures at once. Is passed to :meth:`~Results.get_output_chunk`, defaults to 0.
        """
        block = self.get_output_chunk(begin="{} ATOMIC CHARGES".format(method.upper()), end="-" * 25, match=match)
        ret = []
        for line in block:
            line = line.strip().split()
            if len(line) == 1:  # new set of charges
                ret.append([])
                continue
            if ("Sum" in line) or (not line):  # empty and sum lines
                continue
            ret[-1].append(float(line[-1]))
        if len(ret) == 1:
            return ret[0]
        else:
            return ret

    def get_hirshfeld(self, return_spin=False, match=0, skip=5):
        """Get Hirshfeld Analysis from Output

        :parameter return_spin: Return a tuple of (charge, spin) instead of just charge.
        :parameter match: Select occurence in the output to use. E.g. when running multiple structures at once. Is passed to :meth:`~Results.get_output_chunk`, defaults to 0.
        :parameter skip: Number of lines after the keyword in the outputfile to be skipped. Don't touch if you don't have trouble with your ORCA versions output.
        """
        block = self.get_output_chunk(begin="HIRSHFELD ANALYSIS", end="TOTAL", match=match)
        ret = []
        j = 0
        for i in range(len(block)):
            line = block[j].strip().split()
            if len(line) == 1:  # new set of charges
                ret.append([])
                j += skip
            elif line:
                if return_spin:
                    ret[-1].append(tuple([float(x) for x in line[-2:]]))
                else:
                    ret[-1].append(float(line[-2]))
            j += 1
            if j >= len(block):
                break
        if len(ret) == 1:
            return ret[0]
        else:
            return ret

    def get_orbital_energies(self, unit="a.u.", return_occupancy=False, match=0):
        """Returns Orbital Energies.

        :parameter return_occupancy: Set to *True* to recieve a tuple (Energy, Occupation) for each MO.
        :parameter match: Select occurence in the output to use. E.g. when running multiple structures at once. Is passed to :meth:`~Results.get_output_chunk`, defaults to 0.
        """
        conv = Units.conversion_ratio("a.u.", unit)
        block = self.get_output_chunk(begin="  NO   OCC          E(Eh)            E(eV) ", end="--", match=match)
        ret = []
        fastFWD = False
        for line in block:
            line = line.strip().split()
            if (len(line) == 4) and (line[0] == "0"):  # new set of energies
                ret.append([])
                fastFWD = False
            elif not len(line) == 4:
                fastFWD = True
                continue
            elif fastFWD:
                continue

            if return_occupancy:
                ret[-1].append((float(line[-2]) * conv, float(line[-3])))
            else:
                ret[-1].append(float(line[-2]) * conv)
        if len(ret) == 1:
            return ret[0]
        else:
            return ret


class ORCAJob(SingleJob):
    """
    A class representing a single computational job with
    `ORCA <https://orcaforum.cec.mpg.de>`_.

    In addition to the arguments of |SingleJob|, |ORCAJob| takes a
    ``copy_files`` argument. ``copy_files`` can be a list or string,
    containing paths to files to be copied to the jobs directory.
    This might e.g. be a molecule, restart files etc. By setting
    ``copy_symlink``, the files are not copied, but symlinked with
    relative links. The same things can be passed using the ``settings``
    instance of the job, i.e. ``self.settings.copy_files`` and
    ``self.settings.copy_symlink``. The former overwrites the latter.
    """

    _result_type = ORCAResults

    def __init__(self, copy_files=None, copy_symlink=False, **kwargs):
        SingleJob.__init__(self, **kwargs)
        if copy_files:
            self.settings.copy_files = copy_files
        if copy_symlink:
            self.settings.copy_symlink = copy_symlink

    def _get_ready(self):
        """Copy files to execution dir if self.copy_files is set."""
        SingleJob._get_ready(self)
        if "copy_files" in self.settings:
            if not isinstance(self.settings.copy_files, list):
                copy_files = [self.settings.copy_files]
            else:
                copy_files = self.settings.copy_files
            for f in copy_files:
                if ("copy_symlink" in self.settings) and (self.settings.copy_symlink):
                    symlink(relpath(f, self.path), opj(self.path, basename(f)))
                else:
                    shutil.copy(f, self.path)
        return

    def get_input(self):
        """
        Transform all contents of ``input`` branch of  ``settings`` into
        string with blocks, subblocks, keys and values.
        """

        def get_end(s):
            if (not isinstance(s, Settings)) or ("_end" not in s):
                return s
            else:
                return "{} end".format(s["_end"])

        def pretty_print_inner(s, indent):
            inp = ""
            for i, (key, value) in enumerate(s.items()):
                end = get_end(value)
                if i == 0:
                    inp += " {} {}\n".format(key, end)
                else:
                    inp += "{}{} {}\n".format(indent, key, end)
            return inp

        def pretty_print_orca(s, indent="", print_main=False):
            """Set print_main to true for initial call to have main section
            at the top
            """
            inp = ""
            if print_main:
                inp += "! {}\n\n".format(pretty_print_orca(s.main))
                pretty_print_orca(s, indent)
            if isinstance(s, Settings):
                for k, v in s.items():
                    if k in ("main", "molecule"):
                        # skip the molecule and main section
                        continue
                    else:
                        indent2 = (len(k) + 2) * " "
                        if not isinstance(v, Settings):
                            inp += "%{} {}\n\n".format(k, pretty_print_orca(v))
                        else:
                            block = pretty_print_inner(v, indent2)
                            inp += "%{}{}{}end\n\n".format(k, block, indent)
            elif isinstance(s, list):
                inp += "{}{}".format(indent, " ".join(s))
            else:
                inp += "{}{}".format(indent, s)
            return inp

        inp = pretty_print_orca(self.settings.input, print_main=True)
        if "molecule" in self.settings.input:
            inp += "* {}\n".format(self.settings.input.molecule)
        elif isinstance(self.molecule, Molecule):
            inp += self.print_molecule()
        else:
            raise JobError("ORCA Job needs a molecule to run.")

        return inp

    def print_molecule(self):
        """Print a molecule in the ORCA format using the xyz notation."""
        mol = self.molecule
        if "charge" in mol.properties and isinstance(mol.properties.charge, int):
            charge = mol.properties.charge
        else:
            charge = 0
        if "multiplicity" in mol.properties and isinstance(mol.properties.multiplicity, int):
            multi = mol.properties.multiplicity
        else:
            multi = 1
        # very accurate but this is equal to the accuracy of the ORCA output
        xyz = "\n".join(at.str(symbol=True, space=21, decimal=14) for at in mol.atoms)
        return "* xyz {} {}\n{}\n*\n\n".format(charge, multi, xyz)

    def get_runscript(self):
        """Returned runscript is just one line:
        ``orca myinput.inp``
        """
        return "orca {}".format(self._filename("inp"))

    def check(self):
        """Look for the normal termination signal in the output."""
        s = self.results.grep_output("ORCA TERMINATED NORMALLY")
        return len(s) > 0
