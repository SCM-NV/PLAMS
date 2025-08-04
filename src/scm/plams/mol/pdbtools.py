import numpy as np
from scm.plams.core.errors import FileError, PlamsError
from scm.plams.tools.periodic_table import PeriodicTable

__all__ = ["PDBRecord", "PDBHandler", "PDBAtom"]


_multiline = set(
    [
        "AUTHOR",
        "CAVEAT",
        "COMPND",
        "EXPDTA",
        "MDLTYP",
        "KEYWDS",
        "SOURCE",
        "SPLIT ",
        "SPRSDE",
        "TITLE ",
        "FORMUL",
        "HETNAM",
        "HETSYN",
        "SEQRES",
        "SITE  ",
        "REMARK",
    ]
)

_sequence = [
    "HEADER",
    "OBSLTE",
    "TITLE ",
    "SPLIT ",
    "CAVEAT",
    "COMPND",
    "SOURCE",
    "KEYWDS",
    "EXPDTA",
    "NUMMDL",
    "MDLTYP",
    "AUTHOR",
    "REVDAT",
    "SPRSDE",
    "JRNL  ",
    "REMARK",
    "DBREF ",
    "DBREF1",
    "DBREF2",
    "SEQADV",
    "SEQRES",
    "MODRES",
    "HET   ",
    "HETNAM",
    "HETSYN",
    "FORMUL",
    "HELIX ",
    "SHEET ",
    "SSBOND",
    "LINK  ",
    "CISPEP",
    "SITE  ",
    "CRYST1",
    "ORIGX1",
    "ORIGX2",
    "ORIGX3",
    "SCALE1",
    "SCALE2",
    "SCALE3",
    "MTRIX1",
    "MTRIX2",
    "MTRIX3",
    "MODEL ",
    "CONECT",
    "MASTER",
    "END   ",
]

_coord = ["ATOM  ", "ANISOU", "HETATM", "TER   ", "ENDMDL"]


# ===========================================================================


class PDBRecord:
    __slots__ = ["name", "value", "model"]

    def __init__(self, s=None):
        """
        Instantiate a single record from a PDB file

        * ``s`` -- Either a string representing a line in a PDB file, or one of (PDBAtom, PDBLattice, PDBConnections)
        """
        self.name = None
        self.value = None
        self.model = []

        if s is None:
            return

        if isinstance(s, str):
            s = s.rstrip("\n")
            if len(s) < 80:
                s = "%-80s" % s
            self.name = s[:6]

            if self.name in _coord[:3]:
                self.value = [PDBAtom(s)]
            elif self.name == "CRYST1":
                self.value = [PDBLattice(s)]
            elif self.name == "CONECT":
                self.value = [PDBConnections(s)]
            else:
                self.value = [s[6:]]

        elif isinstance(s, PDBAtom):
            self.name = "ATOM  "
            self.value = [s]

        elif isinstance(s, PDBLattice):
            self.name = "CRYST1"
            self.value = [s]

        elif isinstance(s, PDBConnections):
            self.name = "CONECT"
            self.value = [s]

        else:
            raise PlamsError("No PDB record info supplied")
        self.model = []

    def __str__(self):
        res = []
        if self.name != "_model":
            for val in self.value:
                res += [self.name, str(val), "\n"]
        for i in self.model:
            res += [str(i)]
        res = "".join(res)
        return res

    def is_multiline(self):
        return self.name in _multiline

    def extend(self, s):
        s = s.rstrip("\n")

        def _tonum(ss):
            if ss.isspace():
                return 1
            else:
                return int(ss)

        if self.is_multiline() and self.name == s[:6]:
            val = s[6:]
            if self.name == "REMARK":
                self.value.append(val)
                return True
            beg, end = 1, 4
            if self.name == "FORMUL":
                beg, end = 10, 12
            last = _tonum(self.value[-1][beg:end])
            new = _tonum(val[beg:end])
            if new == last + 1:
                self.value.append(val)
                return True
        return False

    def set_index(self, ind):
        """
        Set the atom index for this record, if relevant
        """
        for val in self.value:
            if isinstance(val, PDBAtom):
                val.index = ind

    def has_atom(self):
        """
        Whether this record represents an atom
        """
        return self.name in ["ATOM  ", "HETATM"]

    def get_atom(self):
        """
        If this record represents an atom, return the PDBAtom object
        """
        if not self.has_atom():
            return None
        return self.value[0]


# ===========================================================================


class PDBHandler:
    def __init__(self, textfile=None):
        self.records = {}
        for key in _sequence + _coord:
            self.records[key] = []
        if textfile is not None:
            if isinstance(textfile, str):
                try:
                    f = open(textfile, "r")
                except:
                    raise FileError("PDBHandler: Error reading file %s" % textfile)
                self.read(f)
                f.close()
            else:  # textfile is an open file object
                self.read(textfile)

    def singlemodel(self):
        return "_model" in self.records

    def read(self, f):
        model = None
        line = f.readline()
        while line:
            newrecord = PDBRecord(line)
            line = f.readline()
            while newrecord.extend(line):
                line = f.readline()
            key = newrecord.name
            self.records[key].append(newrecord)

            if model is not None and key in _coord:
                newrecord.set_index(len(model))
                model.append(newrecord)

            if key == "MODEL ":
                model = newrecord.model
            if model is None and key in _coord[:3]:
                newrecord.set_index(0)
                tmp = PDBRecord("_model")
                model = tmp.model
                model.append(newrecord)
                self.records["_model"] = tmp
            elif key == "END   ":
                break

    def write(self, f):
        for key in _sequence:
            if key == "MODEL " and self.singlemodel():
                f.write(str(self.records["_model"]))
            else:
                for record in self.records[key]:
                    f.write(str(record))

    def calc_master(self):
        def total(key):
            if key in _multiline:
                return sum([len(x.value) for x in self.records[key]])
            else:
                return len(self.records[key])

        remark = total("REMARK")
        het = total("HET   ")
        helix = total("HELIX ")
        sheet = total("SHEET ")
        site = total("SITE  ")
        conect = total("CONECT")
        seqres = total("SEQRES")
        xform = sum(
            map(total, ["ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3"])
        )

        if self.singlemodel():
            ter = total("TER   ")
            coord = total("ATOM  ") + total("HETATM")
        else:
            lst = self.records["MODEL "][0].model
            ter, coord = 0, 0
            for i in lst:
                if i.name == "TER   ":
                    ter += 1
                elif i.name in ["ATOM  ", "HETATM"]:
                    coord += 1

        master = "MASTER    %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i          \n" % (
            remark,
            0,
            het,
            helix,
            sheet,
            0,
            site,
            xform,
            coord,
            ter,
            conect,
            seqres,
        )
        return PDBRecord(master)

    def check_master(self):
        if self.records["MASTER"]:
            old = self.records["MASTER"][0]
            new = self.calc_master()
            return old.value == new.value
        return False

    def get_models(self):
        if self.singlemodel():
            return [self.records["_model"].model]
        else:
            return [x.model for x in self.records["MODEL "]]

    def get_atoms(self, imodel=1):
        """
        Get a list of atom objects for the defined model
        """
        models = self.get_models()
        if imodel > len(models):
            raise PlamsError("There are only %i geometries in the PDB" % (len(models)))

        atoms = []
        for record in models[imodel - 1]:
            if record.has_atom():
                atoms.append(record.get_atom())
        return atoms

    def get_lattice(self):
        """
        Get the lattice as a set of one, two, or three, lattice vectors
        """
        if len(self.records["CRYST1"]) == 0:
            return []
        vectors = self.records["CRYST1"][0].value[0].get_vectors()
        if (vectors[2] ** 2).sum() < 1e-10:
            vectors = vectors[:2]
            if (vectors[1] ** 2).sum() < 1e-10:
                vectors = vectors[:1]
        vectors = vectors.tolist()
        return vectors

    def get_connections(self):
        """
        Get the full connection table
        """
        connections = {}
        for record in self.records["CONECT"]:
            conect = record.value[0]
            if conect.atom_index not in connections:
                connections[conect.atom_index] = []
            for iat in conect.neighbors:
                if iat not in connections[conect.atom_index]:
                    connections[conect.atom_index].append(iat)
        return connections

    def add_record(self, record):
        if record.name in self.records:
            self.records[record.name].append(record)
        elif record.name == "_model":
            self.records[record.name] = record
        else:
            raise PlamsError("PDBHandler.add_record: Invalid record passed: %s" % (record.name))

    def add_model(self, model):
        """
        model: list of PDBRecords of type in _coord
        """
        # Make sure atoms are properly numbered
        for i, record in enumerate(model):
            record.set_index(i)

        # Add the model
        if "_model" in self.records:
            old = self.records["_model"]
            del self.records["_model"]

            newmodel = PDBRecord("MODEL     %4i" % 1)
            newmodel.model = old.model
            endmdl = PDBRecord("ENDMDL")
            newmodel.model.append(endmdl)
            self.add_record(newmodel)
            self.add_record(endmdl)

        if self.records["MODEL "]:  # there were 1+ models present before
            newmodel = PDBRecord("MODEL     %4i" % (1 + len(self.records["MODEL "])))
            newmodel.model = model
            if newmodel.model[-1].name != "ENDMDL":
                newmodel.model.append(PDBRecord("ENDMDL"))
            self.add_record(newmodel)
            self.records["NUMMDL"] = [PDBRecord("NUMMDL    %4i" % len(self.records["MODEL "]))]
        else:
            newmodel = PDBRecord("_model")
            newmodel.model = model
            self.add_record(newmodel)

        for rec in model:
            self.add_record(rec)

        # Make sure that the basic lines are present
        if len(self.records["HEADER"]) == 0:
            self.add_record(PDBRecord("HEADER"))
        if len(self.records["END   "]) == 0:
            self.add_record(PDBRecord("END   "))
        if len(self.records["MASTER"]) == 0:
            self.add_record(self.calc_master())

    def add_atom(self, atom, imodel=1):
        """
        Add an atom to the specified model

        * ``atom`` -- PDBAtom instance
        """
        models = self.get_models()
        nmodels = len(models)
        if nmodels < imodel - 1:
            raise PlamsError("Model %i does not exist" % (imodel))

        record = PDBRecord(atom)
        if nmodels == imodel - 1:
            model = [record]
            self.add_model(model)
        else:
            record.set_index(len(models[imodel - 1]))
            if self.singlemodel():
                self.records["_model"].model.append(record)
            else:
                self.records["MODEL "][imodel - 1].model.append(record)
            self.add_record(record)

    def set_lattice(self, vectors):
        """
        Set the lattice info

        * ``vectors`` -- List of one, two or three vectors
        """
        lattice = PDBLattice.from_vectors(vectors)
        record = PDBRecord(lattice)
        if "CRYTS1" in self.records.keys():
            if len(self.records["CRYST1"]) > 0:
                self.records["CRYST1"] = []
        self.add_record(record)

    def set_connections(self, connection_table):
        """
        Add new connections to the PDB object

        * ``connection_table`` -- Dictionary with atom index as key, and list of neighboring atom indices as value
                                  All indexing starts at 0
        """
        for iat, neighbors in connection_table.items():
            conect = PDBConnections()
            conect.atom_index = iat
            conect.neighbors = neighbors
            record = PDBRecord(conect)
            self.add_record(record)


class PDBAtom:
    """
    Class representing a PDB atom
    """

    def __init__(self, line=None):
        """
        Instantiates an instance of the PDBAtom object
        """
        self.name = "    "
        self.index = 0
        self.coords = None
        self.res = "    "
        self.resnum = "    "
        self.occ = "    "
        self.fix = "    "
        self.seg = ""
        self.element = ""

        if line is not None:
            self._readline(line)

    def _readline(self, line):
        """
        Turn string into atom object
        """
        element = ""

        words = line.split()
        if words[0] not in ["ATOM", "HETATM"]:
            raise PlamsError("This line does not represent an atom")

        # Correct for empty atom name + residue info (as only PLAMS seems to write it)
        if len(line[11:30].strip()) == 0:
            words = words[:2] + ["    ", "    ", "    "] + words[2:]
        # Correct for empty fix and obb values (as only PLAMS seems to write it)
        empty_fixocc = False
        if len(line[54:66].strip()) == 0 and len(words) >= 9:
            nlast = len(line[67:].split())
            words = words[:-nlast] + ["1.00", "0.00"] + words[-nlast:]
            empty_fixocc = True

        # Correct for lack of space between atom name and resname
        if line[13:19] in words[2]:
            words = words[:3] + [line[16:21]] + words[3:]
        words[2] = line[12:16]
        # If there are a segname and en element name, remove and store the latter
        if len(words) > 11:
            if not words[-2].replace(".", "").isdigit():
                element = words[-1]
                words = words[:-1]
        # Take care of spaces in resnum
        if len(words) > 11:
            words[4] = line[21:26]
            words = words[:5] + words[6:]
        # Correct for lack of space between resname and resnum
        if len(words[3]) > 4:
            if len(words[3]) < 9:
                words[3] = words[3][:4]
            else:
                words = words[:3] + [line[17:21], line[22:26]] + words[4:]
        # An extra X after the resnum?
        if words[4] == "X":
            words = words[:3] + [line[17:21], line[22:26]] + words[6:]
        # Find coordinates
        are_coords = [w.replace(".", "") for w in words[5:8]]
        are_coords = [w.replace("-", "") for w in are_coords]
        are_coords = [w.isnumeric() for w in are_coords]
        if False not in are_coords:
            x = float(words[5])
            y = float(words[6])
            z = float(words[7])
        else:
            seg = words[-1]
            b = words[-2]
            o = words[-3]
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            words = words[:5] + [x, y, z, o, b, seg]
        # If there is no segment, add an empty one
        if len(words) == 10:
            words.append("")

        # Check that fix and occ are ok
        if empty_fixocc:
            words[8] = "    "
            words[9] = "    "
        if not words[8].replace(".", "").isnumeric():
            if len(words[8].strip()) > 0:
                raise PlamsError("Bad pdb format\n%s" % (line))
        if not words[9].replace(".", "").isnumeric():
            if len(words[9].strip()) > 0:
                print(words[9])
                raise PlamsError("Bad pdb format\n%s" % (line))

        # Start assigning instance variables
        x = float(words[5])
        y = float(words[6])
        z = float(words[7])
        self.coords = (x, y, z)
        self.name = words[2]
        self.res = words[3]
        self.resnum = words[4]
        self.occ = words[8]
        self.fix = words[9]
        self.seg = words[10]
        self.element = element
        if self.iselement(self.seg):
            self.element = self.seg
            self.seg = ""

    def __str__(self):
        """
        Return as string
        """
        # Define some formatting detauls
        lengths = [5, 4, 4, 6, 8, 8, 8, 6, 6]
        form = ["%5s", " %4s", "%4s", "%6s    ", "%8.3f", "%8.3f", "%8.3f", "%6s", "%6s"]
        form2 = ["%5s", " %4s", " %4s", "%5s    ", "%8.3f", "%8.3f", "%8.3f", "%6s", "%6s"]
        len_seg = 8

        # Check if the residue number is reasonable
        resnum = self.resnum
        if format == "molden":
            if not resnum.replace("-", "").isdigit() and len(resnum.strip()) > 0:
                if not resnum[1:].isdigit():
                    resnum = "0"

        # Prepare the text
        crds = self.coords
        line = [form[1] % ("%i" % (self.index + 1)), self.name, self.res, resnum]
        line += [crds[0], crds[1], crds[2], self.occ, self.fix, self.seg]

        # Check that the coordinates are not too large
        problem = False
        for x in crds:
            if x >= 10000.0 or x <= -1000.0:
                problem = True
        if problem:
            print("Error: Coordinate values too large! ", self.name)

        # Assemble the whole line, except the segment info
        block = ""
        if len(self.res) < 4:
            for f, length, word in zip(form, lengths, line):
                tmpword = word
                if isinstance(word, str):
                    tmpword = word[:length]
                text = f % (tmpword)
                block += text
        else:
            for f, length, word in zip(form2, lengths, line):
                tmpword = word
                if isinstance(word, str):
                    tmpword = word[:length]
                text = f % (tmpword)
                block += text

        # Add the segment and element info
        if len(self.seg.strip()) > 0 or len(self.element.strip()) > 0:
            tmpword = self.seg
            tmpword = tmpword[:len_seg]
            text = "%10s" % (tmpword)
            if len(self.element.strip()) == 2:
                text = "%9s" % (tmpword)
            block += text
            if len(self.element.strip()) > 0:
                text = " %-2s " % (self.element)
                block += text

        return block

    @staticmethod
    def iselement(word):
        """
        Check if word represents element
        """
        elements = [l[0] for l in PeriodicTable.data]
        if word.strip().capitalize() in elements:
            return True
        return False

    def get_symbol(self):
        """
        Try to deduce atomic symbol from e.g. an atom name
        """
        elements = [l[0] for l in PeriodicTable.data]

        word = self.element
        if len(word) == 0:
            word = self.name
        s = word.strip()[:2].capitalize()
        if s[0] in elements:
            s = s[0]
        else:
            if not s in elements:
                s = None
        return s


class PDBLattice:
    """
    Class representing the preiodic lattice of the PDB system
    """

    def __init__(self, line=None):
        """
        Instantiates a PDBLattice object

        * ``line`` -- String representing the PDB line starting with CRYST1
        """
        self.lengths = []
        self.angles = []

        if line is not None:
            values = [float(x) for x in line[6:54].split()]
            self.lengths = values[:3]
            self.angles = values[3:]

    def __str__(self):
        """
        The object as a string
        """
        parts = ["%9.3f%9.3f%9.3f" % tuple(self.lengths)]
        parts += ["%7.2f%7.2f%7.2f P 1           1" % tuple(self.angles)]
        return "".join(parts)

    def get_vectors(self):
        """
        Returns lattice as cell vectors (a 3x3 Numpy array)
        """
        a = self.lengths[0]
        b = self.lengths[1]
        c = self.lengths[2]
        alpha = self.angles[0] * np.pi / 180.0
        beta = self.angles[1] * np.pi / 180.0
        gamma = self.angles[1] * np.pi / 180.0

        va = [a, 0.0, 0.0]
        vb = [b * np.cos(gamma), b * np.sin(gamma), 0.0]

        cx = c * np.cos(beta)
        cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) * c / np.sin(gamma)
        volume = 1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2
        volume += 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
        volume = np.sqrt(volume)
        cz = c * volume / np.sin(gamma)
        vc = [cx, cy, cz]

        cv = np.array([va, vb, vc])
        return cv

    @classmethod
    def from_vectors(cls, vectors):
        """
        Create lattice object from vectors
        """

        def get_angle(va, vb):
            """
            Return angle between two numpy vectors
            """
            lena = np.sqrt((va**2).sum())
            lenb = np.sqrt((vb**2).sum())
            if lenb < 1e-10:
                return 90.0
            cosphi = (va * vb).sum() / (lena * lenb)
            if cosphi > 1:
                cosphi = 1.0
            phi = 360 * (np.arccos(cosphi) / (2 * np.pi))
            return phi

        ret = cls()

        nvecs = len(vectors)
        vecs = np.zeros(9).reshape((3, 3))
        vecs[:nvecs] = vectors

        ret.lengths = np.sqrt((vecs**2).sum(axis=1))
        if ret.lengths.sum() == 0:
            ret.angles = [90.0, 90.0, 90.0]

        alpha = get_angle(vecs[1], vecs[2])
        beta = get_angle(vecs[0], vecs[2])
        gamma = get_angle(vecs[0], vecs[1])
        ret.angles = [alpha, beta, gamma]
        return ret


class PDBConnections:
    """
    Class representing the connecivity of a single atom
    """

    def __init__(self, line=None):
        """
        Instantiates an instance of the PDBConnections class

        * ``line`` -- String from a PDB file, starting with CONECT
        """
        self.atom_index = None
        self.neighbors = []

        if line is None:
            return

        # Extract the connection data
        words = line.split()

        name = words[0]
        if len(name) == 6:
            atnum = int(words[1])
            condata = [int(w) for w in words[2:]]
        elif name[:6] == "CONECT":
            # This suggests that there are no spaces between the atom indices
            condata = []
            factor = 5
            if len(name[6:]) % 6 == 0:
                factor = 6
            atnum = int(name[6 : 6 + factor])
            for i in range((len(name) - 6 - 6) / factor):
                num = int(name[6 + 6 + (factor * i) : 6 + (factor * (i + 1))])
                condata.append(num)

        self.atom_index = atnum - 1
        self.neighbors = [i - 1 for i in condata]

    def __str__(self):
        """
        Write PDBConnections instance as string
        """
        words = [""]
        words += ["%4i" % (self.atom_index + 1)]
        words += ["%4i" % (bat + 1) for bat in self.neighbors]
        return " ".join(words)
