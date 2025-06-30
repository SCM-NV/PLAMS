from typing import List, Literal, Optional, overload, TYPE_CHECKING, Sequence, Dict, Any
import random
import sys
import copy
from warnings import warn

from scm.plams.core.functions import add_to_class, log, requires_optional_package
from scm.plams.mol.atom import Atom
from scm.plams.mol.bond import Bond
from scm.plams.mol.molecule import Molecule
from scm.plams.core.errors import PlamsError

if TYPE_CHECKING:
    from rdkit import Mol as RDKitMol

__all__ = [
    "add_Hs",
    "apply_reaction_smarts",
    "apply_template",
    "gen_coords_rdmol",
    "get_backbone_atoms",
    "modify_atom",
    "to_rdmol",
    "from_rdmol",
    "from_sequence",
    "from_smiles",
    "from_smarts",
    "to_smiles",
    "partition_protein",
    "readpdb",
    "writepdb",
    "get_substructure",
    "get_conformations",
    "yield_coords",
    "canonicalize_mol",
    "to_image",
    "get_reaction_image",
]


@requires_optional_package("rdkit")
def from_rdmol(rdkit_mol: "RDKitMol", confid: int = -1, properties: bool = True) -> Molecule:
    """
    Translate an RDKit molecule into a PLAMS molecule type.
    RDKit properties will be unpickled if their name ends with '_pickled'.

    :parameter rdkit_mol: RDKit molecule
    :type rdkit_mol: rdkit.Chem.Mol
    :parameter int confid: conformer identifier from which to take coordinates
    :parameter bool properties: If all Chem.Mol, Chem.Atom and Chem.Bond properties should be converted from RDKit to PLAMS format.
    :return: a PLAMS molecule
    :rtype: |Molecule|
    """
    from rdkit import Chem

    if isinstance(rdkit_mol, Molecule):
        return rdkit_mol
    # Create PLAMS molecule
    plams_mol = Molecule()
    total_charge = 0
    try:
        Chem.Kekulize(rdkit_mol)
    except Exception:
        pass
    conf = rdkit_mol.GetConformer(id=confid)

    # Add atoms and assign properties to the PLAMS atom if *properties* = True
    for rd_atom in rdkit_mol.GetAtoms():
        pos = conf.GetAtomPosition(rd_atom.GetIdx())
        ch = rd_atom.GetFormalCharge()
        pl_atom = Atom(rd_atom.GetAtomicNum(), coords=(pos.x, pos.y, pos.z), rdkit={"charge": ch})
        if properties and rd_atom.GetPDBResidueInfo():
            pl_atom.properties.rdkit.pdb_info = get_PDBResidueInfo(rd_atom)
        plams_mol.add_atom(pl_atom)
        total_charge += ch

        # Check for R/S information
        stereo = str(rd_atom.GetChiralTag())
        if stereo == "CHI_TETRAHEDRAL_CCW":
            pl_atom.properties.rdkit.stereo = "counter-clockwise"
        elif stereo == "CHI_TETRAHEDRAL_CW":
            pl_atom.properties.rdkit.stereo = "clockwise"

    # Add bonds to the PLAMS molecule
    for bond in rdkit_mol.GetBonds():
        at1 = plams_mol.atoms[bond.GetBeginAtomIdx()]
        at2 = plams_mol.atoms[bond.GetEndAtomIdx()]
        plams_mol.add_bond(Bond(at1, at2, bond.GetBondTypeAsDouble()))

        # Check for cis/trans information
        stereo, bond_dir = str(bond.GetStereo()), str(bond.GetBondDir())
        if stereo == "STEREOZ" or stereo == "STEREOCIS":
            plams_mol.bonds[-1].properties.rdkit.stereo = "Z"
        elif stereo == "STEREOE" or stereo == "STEREOTRANS":
            plams_mol.bonds[-1].properties.rdkit.stereo = "E"
        elif bond_dir == "ENDUPRIGHT":
            plams_mol.bonds[-1].properties.rdkit.stereo = "up"
        elif bond_dir == "ENDDOWNRIGHT":
            plams_mol.bonds[-1].properties.rdkit.stereo = "down"

    # Set charge and assign properties to PLAMS molecule and bonds if *properties* = True
    plams_mol.properties.charge = total_charge
    if properties:
        prop_from_rdmol(plams_mol, rdkit_mol)
        for rd_atom, plams_atom in zip(rdkit_mol.GetAtoms(), plams_mol):
            prop_from_rdmol(plams_atom, rd_atom)
        for rd_bond, plams_bond in zip(rdkit_mol.GetBonds(), plams_mol.bonds):
            prop_from_rdmol(plams_bond, rd_bond)
    return plams_mol


@requires_optional_package("rdkit")
def to_rdmol(
    plams_mol: Molecule,
    sanitize: bool = True,
    properties: bool = True,
    assignChirality: bool = False,
    presanitize: bool = False,
) -> "RDKitMol":
    """
    Translate a PLAMS molecule into an RDKit molecule type.
    PLAMS |Molecule|, |Atom| or |Bond| properties are pickled if they are neither booleans, floats,
    integers, floats nor strings, the resulting property names are appended with '_pickled'.

    :parameter plams_mol: A PLAMS molecule
    :parameter bool sanitize: Kekulize, check valencies, set aromaticity, conjugation and hybridization
    :parameter bool properties: If all |Molecule|, |Atom| and |Bond| properties should be converted from PLAMS to RDKit format.
    :parameter bool assignChirality: Assign R/S and cis/trans information, insofar as this was not yet present in the PLAMS molecule.
    :parameter bool presanitize: Iteratively adjust bonding and atomic charges, to avoid failure of sanitization.
                                 Only relevant is sanitize is set to True.
    :type plams_mol: |Molecule|
    :return: an RDKit molecule
    :rtype: rdkit.Chem.Mol
    """
    from rdkit import Chem, Geometry

    if isinstance(plams_mol, Chem.Mol):
        return plams_mol
    # Create rdkit molecule
    e = Chem.EditableMol(Chem.Mol())

    # Add atoms and assign properties to the RDKit atom if *properties* = True
    for pl_atom in plams_mol.atoms:
        rd_atom = Chem.Atom(int(pl_atom.atnum))
        if "rdkit" in pl_atom.properties:
            if "charge" in pl_atom.properties.rdkit:
                rd_atom.SetFormalCharge(pl_atom.properties.rdkit.charge)
        if properties:
            if "rdkit" in pl_atom.properties:
                if "pdb_info" in pl_atom.properties.rdkit:
                    set_PDBresidueInfo(rd_atom, pl_atom.properties.rdkit.pdb_info)
                for prop in pl_atom.properties.rdkit:
                    if prop not in ("charge", "pdb_info", "stereo"):
                        prop_to_rdmol(rd_atom, prop, pl_atom.properties.rdkit.get(prop))
            prop_dic = {}
            for prop in pl_atom.properties:
                if prop != "rdkit":
                    prop_dic[prop] = pl_atom.properties.get(prop)
            if len(prop_dic) > 0:
                prop_to_rdmol(rd_atom, "plams", prop_dic)

        # Check for R/S information
        if pl_atom.properties.rdkit.stereo:
            stereo = pl_atom.properties.rdkit.stereo.lower()
            if stereo == "counter-clockwise":
                rd_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif stereo == "clockwise":
                rd_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
        e.AddAtom(rd_atom)

    # Mapping of PLAMS bond orders to RDKit bond types:
    def plams_to_rd_bonds(bo):
        if 1.4 < bo < 1.6:
            return 12  # bond type for aromatic bond
        else:
            return int(bo)

    # Add bonds to the RDKit molecule
    for bond in plams_mol.bonds:
        a1 = plams_mol.atoms.index(bond.atom1)
        a2 = plams_mol.atoms.index(bond.atom2)
        e.AddBond(a1, a2, Chem.BondType(plams_to_rd_bonds(bond.order)))
    rdmol = e.GetMol()

    # Check for cis/trans information
    for pl_bond, rd_bond in zip(plams_mol.bonds, rdmol.GetBonds()):
        if pl_bond.properties.rdkit.stereo:
            stereo = pl_bond.properties.rdkit.stereo.lower()
            if stereo == "e" or stereo == "trans":
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOE)
            elif stereo == "z" or stereo == "cis":
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)
            elif stereo == "up":
                rd_bond.SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)
            elif stereo == "down":
                rd_bond.SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)

    # Assign properties to RDKit molecule and bonds if *properties* = True
    # All properties will be taken from 'rdkit' subsettings, except the molecular charge
    if properties:
        prop_dic = {}
        for prop in plams_mol.properties:
            if prop == "rdkit":
                for rdprop in plams_mol.properties.rdkit:
                    prop_to_rdmol(rdmol, rdprop, plams_mol.properties.rdkit.get(rdprop))
            else:
                # prop_dic[prop] = {'plams':plams_mol.properties.get(prop)}
                prop_dic[prop] = plams_mol.properties.get(prop)
        if len(prop_dic) > 0:
            prop_to_rdmol(rdmol, "plams", prop_dic)
        prop_dic = {}
        for pl_bond, rd_bond in zip(plams_mol.bonds, rdmol.GetBonds()):
            for prop in pl_bond.properties:
                if prop == "rdkit":
                    for rdprop in pl_bond.properties.rdkit:
                        if rdprop != "stereo":
                            prop_to_rdmol(rd_bond, rdprop, pl_bond.properties.rdkit.get(rdprop))
                else:
                    prop_dic[prop] = pl_bond.properties.get(prop)
        if len(prop_dic) > 0:
            prop_to_rdmol(rd_bond, "plams", prop_dic)

    if sanitize:
        try:
            if presanitize:
                rdmol = _presanitize(plams_mol, rdmol)
            else:
                Chem.SanitizeMol(rdmol)
        except ValueError as exc:
            log("RDKit Sanitization Error.")
            text = "Most likely this is a problem with the assigned bond orders: "
            text += "Use chemical insight to adjust them."
            log(text)
            log("Note that the atom indices below start at zero, while the AMS-GUI indices start at 1.")
            raise exc
    conf = Chem.Conformer()
    for i, atom in enumerate(plams_mol.atoms):
        xyz = Geometry.Point3D(atom.x, atom.y, atom.z)
        conf.SetAtomPosition(i, xyz)
    rdmol.AddConformer(conf)
    # REB: Assign all stereochemistry, if it wasn't already there
    if assignChirality:
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(rdmol, confId=conf.GetId(), replaceExistingTags=False)
        try:
            Chem.AssignStereochemistryFrom3D(rdmol, confId=conf.GetId(), replaceExistingTags=False)
        except AttributeError:
            pass
    return rdmol


@requires_optional_package("rdkit")
def to_smiles(plams_mol: Molecule, short_smiles: bool = True, **kwargs) -> str:
    """
    Returns the RDKit-generated SMILES string of a PLAMS molecule.

    Note: SMILES strings are generated based on the molecule's connectivity. If the input PLAMS molecule does not contain any bonds, "guessed bonds" will be used.

    :parameter plams_mol: A PLAMS |Molecule|
    :parameter bool short_smiles: whether or not to use some RDKit sanitization to get shorter smiles (e.g. for a water molecule, short_smiles=True -> "O", short_smiles=False -> [H]O[H])
    :parameter \**kwargs: With 'kwargs' you can provide extra optional parameters to the rdkit.Chem method 'MolToSmiles'. See the rdkit documentation for more info.

    :return: the SMILES string
    """
    from rdkit import Chem

    if len(plams_mol.bonds) > 0:
        mol_with_bonds = plams_mol
    else:
        mol_with_bonds = plams_mol.copy()
        mol_with_bonds.guess_bonds()

    rd_mol = to_rdmol(mol_with_bonds, sanitize=False)

    # This sanitization black magic is needed for getting the "short, nice and clean" SMILES string.
    # Without this, the SMILES string for water would be "[H]O[H]". With this is just "O"
    if short_smiles:
        s = Chem.rdmolops.SanitizeFlags
        rdkitSanitizeOptions = (
            s.SANITIZE_ADJUSTHS
            or s.SANITIZE_CLEANUP
            or s.SANITIZE_CLEANUPCHIRALITY
            or s.SANITIZE_FINDRADICALS
            or s.SANITIZE_PROPERTIES
            or s.SANITIZE_SETAROMATICITY
            or s.SANITIZE_SETCONJUGATION
            or s.SANITIZE_SETHYBRIDIZATION
            or s.SANITIZE_SYMMRINGS
        )
        Chem.rdmolops.AssignRadicals(rd_mol)
        rd_mol = Chem.rdmolops.RemoveHs(rd_mol, updateExplicitCount=True, sanitize=False)
        Chem.rdmolops.SanitizeMol(rd_mol, rdkitSanitizeOptions)
    smiles = Chem.MolToSmiles(rd_mol, **kwargs)
    return smiles


pdb_residue_info_items = [
    "AltLoc",
    "ChainId",
    "InsertionCode",
    "IsHeteroAtom",
    "Name",
    "Occupancy",
    "ResidueName",
    "ResidueNumber",
    "SecondaryStructure",
    "SegmentNumber",
    "SerialNumber",
    "TempFactor",
]
# 'MonomerType' was excluded because it is an rdkit type that cannot easilty be serialized


def get_PDBResidueInfo(rdkit_atom):
    pdb_info = {}
    for item in pdb_residue_info_items:
        get_function = "Get" + item
        pdb_info[item] = rdkit_atom.GetPDBResidueInfo().__getattribute__(get_function)()
    return pdb_info


@requires_optional_package("rdkit")
def set_PDBresidueInfo(rdkit_atom, pdb_info):
    from rdkit import Chem

    atom_pdb_residue_info = Chem.AtomPDBResidueInfo()
    for item, value in pdb_info.items():
        set_function = "Set" + item
        atom_pdb_residue_info.__getattribute__(set_function)(value)
    rdkit_atom.SetMonomerInfo(atom_pdb_residue_info)


def prop_to_rdmol(rd_obj, propkey, propvalue):
    """
    Convert a single PLAMS property into an RDKit property.

    :paramter pl_obj: A PLAMS object.
    :type pl_obj: |Molecule|, |Atom| or |Bond|.
    :parameter rd_obj: An RDKit object.
    :type rd_obj: rdkit.Chem.Mol, rdkit.Chem.Atom or rdkit.Chem.Bond
    :parameter str propkey: The |Settings| key of the PLAMS property.
    """
    try:
        import dill as pickle
    except ImportError:
        import pickle

    obj = type(propvalue)
    obj_dict = {bool: rd_obj.SetBoolProp, float: rd_obj.SetDoubleProp, int: rd_obj.SetIntProp, str: rd_obj.SetProp}
    if obj_dict.get(obj):
        obj_dict[obj](propkey, propvalue)
    else:
        name = propkey + "_pickled"
        try:
            rd_obj.SetProp(name, pickle.dumps(propvalue, 0).decode())
        except (Exception, pickle.PicklingError):
            pass


def prop_from_rdmol(pl_obj, rd_obj):
    """
    Convert one or more RDKit properties into PLAMS properties.

    :paramter pl_obj: A PLAMS object.
    :type pl_obj: |Molecule|, |Atom| or |Bond|.
    :parameter rd_obj: An RDKit object.
    :type rd_obj: rdkit.Chem.Mol, rdkit.Chem.Atom or rdkit.Chem.Bond
    """
    try:
        import dill as pickle
    except ImportError:
        import pickle

    prop_dict = rd_obj.GetPropsAsDict()
    for propname in prop_dict.keys():
        if propname == "plams_pickled":
            plams_props = pickle.loads(prop_dict[propname].encode())
            if not isinstance(plams_props, dict):
                raise Exception("PLAMS property not properly stored in RDKit")
            for key, value in plams_props.items():
                pl_obj.properties[key] = value
        else:
            if propname == "__computedProps":
                continue
            if "_pickled" not in propname:
                pl_obj.properties.rdkit[propname] = prop_dict[propname]
            else:
                prop = prop_dict[propname]
                propname = propname.rsplit("_pickled", 1)[0]
                propvalue = pickle.loads(prop.encode())
                pl_obj.properties.rdkit[propname] = propvalue


@overload
def from_smiles(
    smiles: str, nconfs: Literal[1] = ..., name: Optional[str] = ..., forcefield: Optional[str] = ..., rms: float = ...
) -> Molecule: ...


@overload
def from_smiles(
    smiles: str, nconfs: int = ..., name: Optional[str] = ..., forcefield: Optional[str] = ..., rms: float = ...
) -> List[Molecule]: ...


@requires_optional_package("rdkit")
def from_smiles(
    smiles: str, nconfs: int = 1, name: Optional[str] = None, forcefield: Optional[str] = None, rms: float = 0.1
):
    """
    Generates PLAMS molecule(s) from a smiles strings.

    :parameter str smiles: A smiles string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization
        and ranking of conformations. The default value None results in skipping of the
        geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for
        removing similar/equivalent conformations
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    from rdkit import Chem

    smiles = str(smiles.split()[0])
    smiles = Chem.CanonSmiles(smiles)
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdkit_mol.SetProp("smiles", smiles)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)


@requires_optional_package("rdkit")
def from_smarts(
    smarts: str, nconfs: int = 1, name: Optional[str] = None, forcefield: Optional[str] = None, rms: float = 0.1
):
    """
    Generates PLAMS molecule(s) from a smarts strings.
    This allows for example to define hydrogens explicitly.
    However it is less suitable for aromatic molecules (use from_smiles in that case).

    :parameter str smarts: A smarts string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    from rdkit import Chem

    smiles = str(smarts.split()[0])
    mol = Chem.MolFromSmarts(smiles)
    Chem.SanitizeMol(mol)
    molecule = Chem.AddHs(mol)
    molecule.SetProp("smiles", smiles)
    return get_conformations(molecule, nconfs, name, forcefield, rms)


@requires_optional_package("rdkit")
def get_conformations(
    mol,
    nconfs=1,
    name=None,
    forcefield=None,
    rms=-1,
    enforceChirality=False,
    useExpTorsionAnglePrefs="default",
    constraint_ats=None,
    EmbedParameters="EmbedParameters",
    randomSeed=1,
    best_rms=-1,
):
    """
    Generates 3D conformation(s) for an rdkit_mol or a PLAMS Molecule

    :parameter mol: RDKit or PLAMS Molecule
    :type mol: rdkit.Chem.Mol or |Molecule|
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :parameter float best_rms: Root Mean Square deviation of best atomic permutation for removing
        similar/equivalent conformations.
    :parameter bool enforceChirality: Enforce the correct chirality if chiral centers are present
    :parameter str useExpTorsionAnglePrefs: Use experimental torsion angles preferences for the conformer generation by rdkit
    :parameter list constraint_ats: List of atom indices to be constrained
    :parameter str EmbedParameters: Name of RDKit EmbedParameters class ('EmbedParameters', 'ETKDG')
    :parameter int randomSeed: The seed for the random number generator. If set to None the generated conformers will be non-deterministic.
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if isinstance(mol, Molecule):
        if not mol.bonds:
            mol.guess_bonds()
        rdkit_mol = to_rdmol(mol, assignChirality=enforceChirality)
    else:
        rdkit_mol = mol

    def constrained_embedding(rdkit_mol, nconfs, param_obj, template_mol, randomSeed):
        """
        Use RDKit ConstrainedEmbed to add conformers to rdkit_mol (EmbedMultipleConfs does not constrain)

        Note: ConstrainedEmbed seems very sensitive to the randomseed supplied
              High seeds lead to badly constrained coordinates
        """
        # ConstrainedEmbed overwrites conformers (does not add)
        # So I need to do that in a separate molecule, and add the resulting confs later.
        rdmol_orig = copy.deepcopy(rdkit_mol)
        if randomSeed:
            seed = randomSeed
            # Will crash with high seed
            if randomSeed > 1e9:
                seed = randomSeed - (int((randomSeed / 1e9)) * int(1e9))
            random_ints = [i + seed for i in range(nconfs)]
        else:
            base = random.getrandbits(32)
            random_ints = [i + base for i in range(nconfs)]
        conformers = []
        for seed in random_ints:
            rdmol_orig.RemoveAllConformers()
            rdmol_orig = AllChem.ConstrainedEmbed(rdmol_orig, template_mol, param_obj, randomseed=seed)
            conf = rdmol_orig.GetConformer(0)
            newconf = Chem.Conformer()
            for i, pos in enumerate(conf.GetPositions()):
                newconf.SetAtomPosition(i, pos)
            conformers.append(newconf)
            rdkit_mol.AddConformer(newconf, True)
            cids = [c.GetId() for c in rdkit_mol.GetConformers()]
        return cids

    def MMFFenergy(cid):
        ff = AllChem.MMFFGetMoleculeForceField(rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol), confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = (
                "MMFF energy calculation failed for molecule: "
                + Chem.MolToSmiles(rdkit_mol)
                + "\nNo geometry optimization was performed."
            )
            warn(msg)
            energy = 1e9
        return energy

    def UFFenergy(cid):
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = (
                "MMFF energy calculation failed for molecule: "
                + Chem.MolToSmiles(rdkit_mol)
                + "\nNo geometry optimization was performed."
            )
            warn(msg)
            energy = 1e9
        return energy

    def remove_some_Hs(m):
        res = Chem.RWMol(m)
        c_hs = [x[0] for x in m.GetSubstructMatches(Chem.MolFromSmarts("[#1;$([#1]-[#6])]"))]
        c_hs.sort(reverse=True)
        for aid in c_hs:
            res.RemoveAtom(aid)
        return res.GetMol()

    if name:
        rdkit_mol.SetProp("name", name)

    if best_rms > 0:
        if rms > 0:
            raise PlamsError("Cannot set both rms and best_rms")
        rms = best_rms

    # if enforceChirality :
    #    # This is how chirality is enforced in the GUI. The argument is not passed to AllChem.EmbedMultipleConfs
    #    Chem.AssignAtomChiralTagsFromStructure(rdkit_mol)
    # param_obj = AllChem.ETKDG()
    param_obj = getattr(AllChem, EmbedParameters)()
    param_obj.pruneRmsThresh = rms
    param_obj.enforceChirality = enforceChirality
    if useExpTorsionAnglePrefs != "default":  # The default (True of False) changes with rdkit versions
        param_obj.useExpTorsionAnglePrefs = True
    if constraint_ats is not None:
        # Just adding the coordMap and using EmbedMultipleConfs does not seem to work (in this version)
        # coordMap = {}
        # for i, iat in enumerate(constraint_ats):
        #     coordMap[iat] = rdkit_mol.GetConformer(0).GetAtomPosition(iat)
        # param_obj.coordMap = coordMap
        emol = Chem.RWMol(rdkit_mol)
        indices = [i for i in range(rdkit_mol.GetNumAtoms()) if not i in constraint_ats][::-1]
        for iat in indices:
            emol.RemoveAtom(iat)
        template_mol = emol.GetMol()
    else:
        param_obj.randomSeed = randomSeed if randomSeed is not None else random.getrandbits(31)
    try:
        if constraint_ats is not None:
            cids = constrained_embedding(rdkit_mol, nconfs, param_obj, template_mol, randomSeed)
        else:
            param_obj.randomSeed = randomSeed if randomSeed is not None else random.getrandbits(31)
            cids = list(AllChem.EmbedMultipleConfs(rdkit_mol, nconfs, param_obj))
    except Exception:
        # ``useRandomCoords = True`` prevents (poorly documented) crash for large systems
        param_obj.useRandomCoords = True
        if constraint_ats is not None:
            cids = constrained_embedding(rdkit_mol, nconfs, param_obj, template_mol, randomSeed)
        else:
            cids = list(AllChem.EmbedMultipleConfs(rdkit_mol, nconfs, param_obj))
    if len(cids) == 0:
        # Sometimes rdkit does not crash (for large systems), but simply doe snot create conformers
        param_obj.useRandomCoords = True
        if constraint_ats is not None:
            cids = constrained_embedding(rdkit_mol, nconfs, param_obj, template_mol, randomSeed)
        else:
            cids = list(AllChem.EmbedMultipleConfs(rdkit_mol, nconfs, param_obj))

    if forcefield:
        # Select the forcefield (UFF or MMFF)
        optimize_molecule, energy = {
            "uff": [AllChem.UFFOptimizeMolecule, UFFenergy],
            "mmff": [AllChem.MMFFOptimizeMolecule, MMFFenergy],
        }[forcefield]

        # Optimize and sort conformations
        for cid in cids:
            optimize_molecule(rdkit_mol, confId=cid)
        cids.sort(key=energy)

    # Remove duplicate conformations based on RMS
    if best_rms > 0 or forcefield:
        rdmol_local = rdkit_mol
        rms_function = AllChem.AlignMol
        if best_rms > 0:
            # Remove the H atoms, and prepare to use the more expensive RDKit function
            rdmol_local = remove_some_Hs(rdkit_mol)
            rms_function = AllChem.GetBestRMS
        keep = [cids[0]]
        for cid in cids[1:]:
            for idx in keep:
                try:
                    # r = AllChem.AlignMol(rdkit_mol, rdkit_mol, cid, idx)
                    r = rms_function(rdmol_local, rdmol_local, cid, idx)
                except Exception:
                    r = rms + 1
                    message = "Alignment failed in multiple conformation generation: "
                    message += Chem.MolToSmiles(rdkit_mol)
                    message += "\nAssuming different conformations."
                    warn(message)
                if r < rms:
                    break
            else:
                keep.append(cid)
        cids = keep

    if nconfs == 1:
        return from_rdmol(rdkit_mol)
    else:
        return [from_rdmol(rdkit_mol, cid) for cid in cids]


@requires_optional_package("rdkit")
def from_sequence(sequence, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates PLAMS molecule from a peptide sequence.
    Includes explicit hydrogens and 3D coordinates.

    :parameter str sequence: A peptide sequence, e.g. 'HAG'
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :return: A peptide molecule with hydrogens and 3D coordinates
        or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    from rdkit import Chem

    rdkit_mol = Chem.AddHs(Chem.MolFromSequence(sequence))
    rdkit_mol.SetProp("sequence", sequence)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)


@requires_optional_package("rdkit")
def calc_rmsd(mol1, mol2):
    """
    Superimpose two molecules and calculate the root-mean-squared deviations of
    the atomic positions.
    The molecules should be identical, but the ordering of the atoms may differ.

    :param mol1: Molecule 1
    :param mol2: Molecule 2
    :return: The rmsd after superposition
    :rtype: float
    """
    from rdkit.Chem import AllChem

    rdkit_mol1 = to_rdmol(mol1)
    rdkit_mol2 = to_rdmol(mol2)
    try:
        return AllChem.GetBestRMS(rdkit_mol1, rdkit_mol2)
    except:
        return -999


@requires_optional_package("rdkit")
def modify_atom(mol, idx, element):
    """
    Change atom "idx" in molecule "mol" to "element" and add or remove hydrogens accordingly

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter int idx: index of the atom to be modified
    :parameter str element:
    :return: Molecule with new element and possibly added or removed hydrogens
    :rtype: |Molecule|
    """
    from rdkit import Chem

    rdmol = to_rdmol(mol)
    if rdmol.GetAtomWithIdx(idx).GetSymbol() == element:
        return mol
    else:
        e = Chem.EditableMol(rdmol)
        for neighbor in reversed(rdmol.GetAtomWithIdx(idx - 1).GetNeighbors()):
            if neighbor.GetSymbol() == "H":
                e.RemoveAtom(neighbor.GetIdx())
        e.ReplaceAtom(idx - 1, Chem.Atom(element))
        newmol = e.GetMol()
        Chem.SanitizeMol(newmol)
        newmol = Chem.AddHs(newmol, addCoords=True)
        return from_rdmol(newmol)


@requires_optional_package("rdkit")
def apply_template(mol, template):
    """
    Modifies bond orders in PLAMS molecule according template smiles structure.

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter str template: smiles string defining the correct chemical structure
    :return: Molecule with correct chemical structure and provided 3D coordinates
    :rtype: |Molecule|
    """
    from rdkit import Chem

    rdmol = to_rdmol(mol, sanitize=False)
    template_mol = Chem.AddHs(Chem.MolFromSmiles(template))
    newmol = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, rdmol)
    return from_rdmol(newmol)


@requires_optional_package("rdkit")
def apply_reaction_smarts(mol, reaction_smarts, complete=False, forcefield=None, return_rdmol=False):
    """
    Applies reaction smirks and returns product.
    If returned as a PLAMS molecule, thismolecule.properties.orig_atoms
    is a list of indices of atoms that have not been changed
    (which can for example be used partially optimize new atoms only with the freeze keyword)

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter str reactions_smarts: Reactions smarts to be applied to molecule
    :parameter complete: Apply reaction until no further changes occur or given
        fraction of reaction centers have been modified
    :type complete: bool or float (value between 0 and 1)
    :parameter forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization of product structures.
    :type forcefield: str
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: (product molecule, list of unchanged atoms)
    :rtype: (|Molecule|, list of int)
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def react(reactant, reaction):
        """Apply reaction to reactant and return products"""
        ps = reaction.RunReactants([reactant])
        # if reaction doesn't apply, return the reactant
        if len(ps) == 0:
            return [(reactant, range(reactant.GetNumAtoms()))]
        full = len(ps)
        while complete:  # when complete is True
            # apply reaction until no further changes
            r = random.randint(0, len(ps) - 1)
            reactant = ps[r][0]
            ps = reaction.RunReactants([reactant])
            if len(ps) == 0 or len(ps) / full < (1 - complete):
                ps = [[reactant]]
                break
        # add hydrogens and generate coordinates for new atoms
        products = []
        for p in ps[0]:
            Chem.SanitizeMol(p)
            q = Chem.AddHs(p)
            Chem.SanitizeMol(q)
            u = gen_coords_rdmol(q)  # These are the atoms that have not changed
            products.append((q, u))
        return products

    mol = to_rdmol(mol)
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    # RDKit removes fragments that are disconnected from the reaction center
    # In order to keep these, the molecule is first split in separate fragments
    # and the results, including non-reacting parts, are re-combined afterwards
    frags = Chem.GetMolFrags(mol, asMols=True)
    product = Chem.Mol()
    unchanged = []  # List of atoms that have not changed
    for frag in frags:
        for p, u in react(frag, reaction):
            unchanged += [product.GetNumAtoms() + i for i in u]
            product = Chem.CombineMols(product, p)
    if forcefield:
        optimize_coordinates(product, forcefield, fixed=unchanged)
    # The molecule is returned together with a list of atom indices of the atoms
    # that are identical to those
    # in the reactants. This list can be used in subsequent partial optimization of the molecule
    if not return_rdmol:
        product = from_rdmol(product)
        product.properties.orig_atoms = [a + 1 for a in unchanged]
    return product


def gen_coords(plamsmol):
    """Calculate 3D positions only for atoms without coordinates"""
    rdmol = to_rdmol(plamsmol)
    unchanged = gen_coords_rdmol(rdmol)
    conf = rdmol.GetConformer()
    for a in range(len(plamsmol.atoms)):
        pos = conf.GetAtomPosition(a)
        atom = plamsmol.atoms[a]
        atom._setx(pos.x)
        atom._sety(pos.y)
        atom._setz(pos.z)
    return [a + 1 for a in unchanged]


@requires_optional_package("rdkit")
def gen_coords_rdmol(rdmol):
    from rdkit.Chem import AllChem

    ref = rdmol.__copy__()
    conf = rdmol.GetConformer()
    coordDict = {}
    unchanged = []
    maps = []
    # Put known coordinates in coordDict
    for i in range(rdmol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if (-0.0001 < pos.x < 0.0001) and (-0.0001 < pos.y < 0.0001) and (-0.0001 < pos.z < 0.0001):
            continue  # atom without coordinates
        coordDict[i] = pos
        unchanged.append(i)
        maps.append((i, i))
    # compute coordinates for new atoms, keeping known coordinates
    rms = 1
    rs = 1
    # repeat embedding and alignment until the rms of mapped atoms is sufficiently small
    if rdmol.GetNumAtoms() > len(maps):
        while rms > 0.1:
            AllChem.EmbedMolecule(rdmol, coordMap=coordDict, randomSeed=rs, useBasicKnowledge=True)
            # align new molecule to original coordinates
            rms = AllChem.AlignMol(rdmol, ref, atomMap=maps)
            rs += 1
    return unchanged


@requires_optional_package("rdkit")
def optimize_coordinates(rdkit_mol, forcefield, fixed=[]):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def MMFFminimize():
        ff = AllChem.MMFFGetMoleculeForceField(rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol))
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("MMFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))

    def UFFminimize():
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, ignoreInterfragInteractions=True)
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("UFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))

    optimize_molecule = {"uff": UFFminimize, "mmff": MMFFminimize}[forcefield]
    Chem.SanitizeMol(rdkit_mol)
    optimize_molecule()
    return


@requires_optional_package("rdkit")
def write_molblock(plams_mol, file=sys.stdout):
    from rdkit import Chem

    file.write(Chem.MolToMolBlock(to_rdmol(plams_mol)))


@requires_optional_package("rdkit")
def readpdb(pdb_file, sanitize=True, removeHs=False, proximityBonding=False, return_rdmol=False):
    """
    Generate a molecule from a PDB file

    :param pdb_file: The PDB file to read
    :type pdb_file: path- or file-like
    :param bool sanitize:
    :param bool removeHs: Hydrogens are removed if True
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: The molecule
    :rtype: |Molecule| or rdkit.Chem.Mol
    """
    from rdkit import Chem

    try:
        pdb_file = open(pdb_file, "r")
    except TypeError:
        pass  # pdb_file is a file-like object... hopefully

    pdb_mol = Chem.MolFromPDBBlock(pdb_file.read(), sanitize=sanitize, removeHs=removeHs)
    return pdb_mol if return_rdmol else from_rdmol(pdb_mol)


@requires_optional_package("rdkit")
def writepdb(mol, pdb_file=sys.stdout):
    """
    Write a PDB file from a molecule

    :parameter mol: molecule to be exported to PDB
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param pdb_file: The PDB file to write to, or a filename
    :type pdb_file: path- or file-like
    """
    from rdkit import Chem

    try:
        pdb_file = open(pdb_file, "w")
    except TypeError:
        pass  # pdb_file is a file-like object... hopefully

    mol = to_rdmol(mol, sanitize=False)
    pdb_file.write(Chem.MolToPDBBlock(mol))


@requires_optional_package("rdkit")
def add_Hs(mol, forcefield=None, return_rdmol=False):
    """
    Add hydrogens to protein molecules read from PDB.
    Makes sure that the hydrogens get the correct PDBResidue info.

    :param mol: Molecule to be protonated
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param str forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization on new atoms.
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: A molecule with explicit hydrogens added
    :rtype: |Molecule| or rdkit.Chem.Mol
    """
    from rdkit import Chem

    mol = to_rdmol(mol)
    retmol = Chem.AddHs(mol)
    for atom in retmol.GetAtoms():
        if atom.GetPDBResidueInfo() is None and atom.GetSymbol() == "H":
            bond = atom.GetBonds()[0]
            if bond.GetBeginAtom().GetIdx() == atom.GetIdx:
                connected_atom = bond.GetEndAtom()
            else:
                connected_atom = bond.GetBeginAtom()
            try:
                ResInfo = connected_atom.GetPDBResidueInfo()
                if ResInfo is None:
                    continue  # Segmentation faults are raised if ResInfo is None
                atom.SetMonomerInfo(ResInfo)
                atomname = "H" + atom.GetPDBResidueInfo().GetName()[1:]
                atom.GetPDBResidueInfo().SetName(atomname)
            except:
                pass
    unchanged = gen_coords_rdmol(retmol)
    if forcefield:
        optimize_coordinates(retmol, forcefield, fixed=unchanged)
    return retmol if return_rdmol else from_rdmol(retmol)


@requires_optional_package("rdkit")
def add_fragment(rwmol, frag, rwmol_atom_idx=None, frag_atom_idx=None, bond_order=None):
    from rdkit import Chem

    molconf = rwmol.GetConformer()
    fragconf = frag.GetConformer()
    new_indices = []
    for a in frag.GetAtoms():
        new_index = rwmol.AddAtom(a)
        new_indices.append(new_index)
        molconf.SetAtomPosition(new_index, fragconf.GetAtomPosition(a.GetIdx()))
    for b in frag.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        rwmol.AddBond(new_indices[ba], new_indices[ea], b.GetBondType())
    if bond_order:
        rwmol.AddBond(rwmol_atom_idx, new_indices[frag_atom_idx], Chem.BondType.values[bond_order])
        rwmol.GetAtomWithIdx(new_indices[frag_atom_idx]).SetNumRadicalElectrons(0)


@requires_optional_package("rdkit")
def get_fragment(mol, indices, incl_expl_Hs=True, neutralize=True):
    from rdkit import Chem

    molconf = mol.GetConformer()
    fragment = Chem.RWMol(Chem.Mol())
    fragconf = Chem.Conformer()
    # Put atoms in fragment
    for i in indices:
        atom = mol.GetAtomWithIdx(i)
        new_index = fragment.AddAtom(atom)
        pos = molconf.GetAtomPosition(i)
        fragconf.SetAtomPosition(new_index, pos)
    # Put bonds in fragment
    for b in mol.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        if ba in indices and ea in indices:
            fragment.AddBond(indices.index(ba), indices.index(ea), b.GetBondType())
            continue
        if not incl_expl_Hs:
            continue
        if ba in indices and mol.GetAtomWithIdx(ea).GetSymbol() == "H":
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ea))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ea))
            fragment.AddBond(indices.index(ba), hi, Chem.BondType.SINGLE)
            continue
        if ea in indices and mol.GetAtomWithIdx(ba).GetSymbol() == "H":
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ba))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ba))
            fragment.AddBond(indices.index(ea), hi, Chem.BondType.SINGLE)
    ret_frag = fragment.GetMol()
    Chem.SanitizeMol(ret_frag)
    if neutralize:
        for atom in ret_frag.GetAtoms():
            nrad = atom.GetNumRadicalElectrons()
            if nrad > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + nrad)
                atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(ret_frag)
    ret_frag.AddConformer(fragconf)
    return ret_frag


@requires_optional_package("rdkit")
def partition_protein(mol, residue_bonds=None, split_heteroatoms=True, return_rdmol=False):
    """
    Splits a protein molecule into capped amino acid fragments and caps.

    :param mol: A protein molecule
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param tuple residue_bonds: a tuple of pairs of residue number indicating which
        peptide bonds to split. If none, split all peptide bonds.
    :param bool split_heteroatoms: if True, all bonds between a heteroatom and
        a non-heteroatom across residues are removed
    :return: list of fragments, list of caps
    """
    from rdkit import Chem

    mol = to_rdmol(mol)
    caps = []
    em = Chem.RWMol(mol)
    if split_heteroatoms:
        for bond in mol.GetBonds():
            resinfa = bond.GetBeginAtom().GetPDBResidueInfo()
            resinfb = bond.GetEndAtom().GetPDBResidueInfo()
            if resinfa.GetIsHeteroAtom() is not resinfb.GetIsHeteroAtom():
                if resinfa.GetResidueNumber() != resinfb.GetResidueNumber():
                    em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    # Split peptide bonds
    pept_bond = Chem.MolFromSmarts("[C;X4;H1,H2][CX3](=O)[NX3][C;X4;H1,H2][CX3](=O)")
    for match in mol.GetSubstructMatches(pept_bond):
        if residue_bonds:
            resa = mol.GetAtomWithIdx(match[1]).GetPDBResidueInfo().GetResidueNumber()
            resb = mol.GetAtomWithIdx(match[3]).GetPDBResidueInfo().GetResidueNumber()
            if (resa, resb) not in residue_bonds and (resb, resa) not in residue_bonds:
                continue
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_o_ind = cap.GetSubstructMatch(Chem.MolFromSmarts("[C;X4][CX3]=O"))
        cap_o = get_fragment(cap, cap_o_ind, neutralize=False)
        cap_n_ind = cap.GetSubstructMatch(Chem.MolFromSmarts("O=[CX3][NX3][C;X4]"))[2:]
        cap_n = get_fragment(cap, cap_n_ind, neutralize=False)
        em.RemoveBond(match[1], match[3])
        add_fragment(em, cap_o, match[3], 1, 1)
        add_fragment(em, cap_n, match[1], 0, 1)
    # Split disulfide bonds
    ss_bond = Chem.MolFromSmarts("[C;X4;H1,H2]SS[C;X4;H1,H2]")
    for match in mol.GetSubstructMatches(ss_bond):
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_s_ind = cap.GetSubstructMatch(Chem.MolFromSmarts("[C;X4]SS[C;X4]"))
        cap_s1 = get_fragment(cap, cap_s_ind[0:2], neutralize=False)
        cap_s2 = get_fragment(cap, cap_s_ind[2:4], neutralize=False)
        em.RemoveBond(match[1], match[2])
        add_fragment(em, cap_s1, match[2], 1, 1)
        add_fragment(em, cap_s2, match[1], 0, 1)
    frags = Chem.GetMolFrags(em.GetMol(), asMols=True, sanitizeFrags=False)
    if not return_rdmol:
        frags = [from_rdmol(frag) for frag in frags]
    return frags, caps


@requires_optional_package("rdkit")
def charge_AAs(mol, return_rdmol=False):
    from rdkit import Chem

    ionizations = {"ARG_NH2": 1, "LYS_NZ": 1, "GLU_OE2": -1, "ASP_OD2": -1}
    mol = to_rdmol(mol)
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        res_atom = resinfo.GetResidueName() + "_" + resinfo.GetName().strip()
        try:
            atom.SetFormalCharge(ionizations[res_atom])
            Chem.SanitizeMol(mol)
        except KeyError:
            pass
        Chem.SanitizeMol(mol)
    return mol if return_rdmol else from_rdmol(mol)


def get_backbone_atoms(mol):
    """
    Return a list of atom indices corresponding to the backbone atoms in a peptide molecule.
    This function assumes PDB information in properties.pdb_info of each atom, which is the case
    if the molecule is generated with the "readpdb" or "from_sequence" functions.

    :parameter mol: a peptide molecule
    :type mol: |Molecule| or rdkit.Chem.Mol
    :return: a list of atom indices
    :rtype: list
    """
    mol = from_rdmol(mol)
    backbone = ["N", "CA", "C", "O"]
    return [a for a in range(1, len(mol) + 1) if str(mol[a].properties.pdb_info.Name).strip() in backbone]


@requires_optional_package("rdkit")
def get_substructure(mol, func_list):
    """
    Search for functional groups within a molecule based on a list of reference functional groups.
    SMILES strings, PLAMS and/or RDKit molecules can be used interchangeably in "func_list".

    Example:

    .. code:: python

        >>> mol = from_smiles('OCCO')  # Ethylene glycol
        >>> func_list = ['[H]O', 'C[N+]', 'O=PO']
        >>> get_substructure(mol, func_list)

        {'[H]O': [(<scm.plams.mol.atom.Atom at 0x125183518>,
                   <scm.plams.mol.atom.Atom at 0x1251836a0>),
                  (<scm.plams.mol.atom.Atom at 0x125183550>,
                   <scm.plams.mol.atom.Atom at 0x125183240>)]}

    :parameter mol: A PLAMS molecule.
    :type mol: |Molecule|
    :parameter list func_list: A list of functional groups.
        Functional groups can be represented by SMILES strings, PLAMS and/or RDKit molecules.
    :return: A dictionary with functional groups from "func_list" as keys and a list of n-tuples
        with matching PLAMS |Atom| as values.
    """
    from rdkit import Chem

    def _to_rdmol(functional_group):
        """Turn a SMILES strings, RDKit or PLAMS molecules into an RDKit molecule."""
        if isinstance(functional_group, str):
            # RDKit tends to remove explicit hydrogens if SANITIZE_ADJUSTHS is enabled
            sanitize = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
            ret = Chem.MolFromSmiles(functional_group, sanitize=False)
            Chem.rdmolops.SanitizeMol(ret, sanitizeOps=sanitize)
            return ret
        elif isinstance(functional_group, Molecule):
            return to_rdmol(functional_group)
        elif isinstance(functional_group, Chem.Mol):
            return functional_group
        raise TypeError(
            "get_substructure: "
            + str(type(functional_group))
            + " is not a supported \
                        object type"
        )

    def _get_match(mol, rdmol, functional_group):
        """Perform a substructure match on "mol".
        If a match is found, return a list of n-tuples consisting PLAMS |Atom|.
        Otherwise return False."""
        matches = rdmol.GetSubstructMatches(functional_group)
        if matches:
            return [tuple(mol[j + 1] for j in idx_tup) for idx_tup in matches]
        return False

    rdmol = to_rdmol(mol)
    rdmol_func_list = [_to_rdmol(i) for i in func_list]
    gen = (_get_match(mol, rdmol, i) for i in rdmol_func_list)
    return {key: value for key, value in zip(func_list, gen) if value}


def yield_coords(rdmol, id=-1):
    """Take an rdkit molecule and yield its coordinates as 3-tuples.

    .. code-block:: python

        >>> from scm.plams import yield_coords
        >>> from rdkit import Chem

        >>> rdmol = Chem.Mol(...)  # e.g. Methane
        >>> for xyz in yield_coords(rdmol):
        ...     print(xyz)
        (-0.0, -0.0, -0.0)
        (0.6405, 0.6405, -0.6405)
        (0.6405, -0.6405, 0.6405)
        (-0.6405, 0.6405, 0.6405)
        (-0.6405, -0.6405, -0.6405)


    The iterator produced by this function can, for example, be passed to
    :meth:`Molecule.from_array()<scm.plams.mol.molecule.Molecule.from_array>`
    the update the coordinates of a PLAMS Molecule in-place.

    .. code-block:: python

        >>> from scm.plams import Molecule

        >>> mol = Molecule(...)

        >>> xyz_iterator = yield_coords(rdmol)
        >>> mol.from_array(xyz_iterator)


    :parameter rdmol: An RDKit mol.
    :type rdmol: rdkit.Chem.Mol
    :parameter int id: The ID of the desired conformer.
    :return: An iterator yielding 3-tuples with *rdmol*'s Cartesian coordinates.
    :rtype: iterator
    """
    conf = rdmol.GetConformer(id=id)
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        yield (pos.x, pos.y, pos.z)


@add_to_class(Molecule)
def assign_chirality(self):
    """
    Assigns stereo-info to PLAMS molecule by invoking RDKIT
    """
    rd_mol = to_rdmol(self, assignChirality=True)
    pl_mol = from_rdmol(rd_mol)

    # Add R/S info to self
    for iat, pl_atom in enumerate(pl_mol.atoms):
        # Check for R/S information
        if pl_atom.properties.rdkit.stereo:
            self.atoms[iat].properties.rdkit.stereo = pl_atom.properties.rdkit.stereo

    # Add cis/trans information to self
    for ibond, pl_bond in enumerate(pl_mol.bonds):
        if pl_bond.properties.rdkit.stereo:
            self.bonds[ibond] = pl_bond.properties.rdkit.stereo


@add_to_class(Molecule)
@requires_optional_package("rdkit")
def get_chirality(self):
    """
    Returns the chirality of the atoms
    """
    from rdkit import Chem

    rd_mol = to_rdmol(self, assignChirality=True)
    return Chem.FindMolChiralCenters(rd_mol, force=True, includeUnassigned=True)


@requires_optional_package("rdkit")
def canonicalize_mol(mol, inplace=False, **kwargs):
    r"""Take a PLAMS molecule and sort its atoms based on their canonical rank.

    Example:

    .. code:: python

        >>> from scm.plams import Molecule, canonicalize_mol

        # Methane
        >>> mol: Molecule = ...
        >>> print(mol)
        Atoms:
            1         H      0.640510      0.640510     -0.640510
            2         H      0.640510     -0.640510      0.640510
            3         C      0.000000      0.000000      0.000000
            4         H     -0.640510      0.640510      0.640510
            5         H     -0.640510     -0.640510     -0.640510

        >>> print(canonicalize_mol(mol))
        Atoms:
            1         C      0.000000      0.000000      0.000000
            2         H     -0.640510     -0.640510     -0.640510
            3         H     -0.640510      0.640510      0.640510
            4         H      0.640510     -0.640510      0.640510
            5         H      0.640510      0.640510     -0.640510

    :parameter mol: The to-be canonicalized molecule.
    :type mol: |Molecule|
    :parameter bool inplace: Whether to sort the atoms inplace or to return a new molecule.
    :parameter \**kwargs: Further keyword arguments for rdkit.Chem.CanonicalRankAtoms_.
    :return: Either ``None`` or a newly sorted molecule, depending on the value of ``inplace``.
    :rtype: None or |Molecule|

    .. _rdkit.Chem.CanonicalRankAtoms: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.CanonicalRankAtoms

    """
    from rdkit import Chem

    if not isinstance(mol, Molecule):
        raise TypeError("`mol` expected a plams Molecule")
    rdmol = to_rdmol(mol)
    idx_rank = Chem.CanonicalRankAtoms(rdmol, **kwargs)

    if inplace:
        mol.atoms = [at for _, at in sorted(zip(idx_rank, mol.atoms), reverse=True)]
        return None
    else:
        ret = mol.copy()
        ret.atoms = [at for _, at in sorted(zip(idx_rank, ret.atoms), reverse=True)]
        return ret


@requires_optional_package("rdkit")
@requires_optional_package("PIL")
def to_image(
    mol: Molecule,
    remove_hydrogens: bool = True,
    filename: Optional[str] = None,
    fmt: str = "svg",
    size: Sequence[int] = (200, 100),
    as_string: bool = True,
):
    """
    Convert single molecule to single image object

    * ``mol`` -- PLAMS Molecule object
    * ``remove_hydrogens`` -- Wether or not to remove the H-atoms from the image
    * ``filename`` -- Optional: Name of image file to be created.
    * ``fmt`` -- One of "svg", "png", "eps", "pdf", "jpeg"
    * ``size`` -- Tuple/list containing width and height of image in pixels.
    * ``as_string`` -- Returns the image as a string or bytestring. If set to False, the original format
                       will be returned, which can be either a PIL image or SVG text
                       We do this because after converting a PIL image to a bytestring it is not possible
                       to further edit it (with our version of PIL).
    * Returns -- Image text file / binary of image text file / PIL Image object.
    """
    from io import BytesIO
    from PIL import Image
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.Draw import MolsToGridImage

    extensions = ["svg", "png", "eps", "pdf", "jpeg"]

    classes: Dict[str, Any] = {}
    classes["svg"] = _MolsToGridSVG
    for ext in extensions[1:]:
        classes[ext] = None
    # PNG can only be created in this way with later version of RDKit
    if hasattr(rdMolDraw2D, "MolDraw2DCairo"):
        for ext in extensions[1:]:
            classes[ext] = MolsToGridImage

    # Determine the type of image file
    if filename is not None:
        if "." in filename:
            extension = filename.split(".")[-1]
            if extension in classes.keys():
                fmt = extension
            else:
                msg = [f"Image type {extension} not available."]
                msg += [f"Available extensions are: {' '.join(extensions)}"]
                raise Exception("\n".join(msg))
        else:
            filename = ".".join([filename, fmt])

    if fmt not in classes.keys():
        raise Exception(f"Image type {fmt} not available.")

    rdmol = _rdmol_for_image(mol, remove_hydrogens)

    # Draw the image
    if classes[fmt.lower()] is None:
        # With AMS version of RDKit MolsToGridImage fails for eps (because of paste)
        img = Draw.MolToImage(rdmol, size=size)
        buf = BytesIO()
        img.save(buf, format=fmt)
        img_text = buf.getvalue()
    else:
        # This fails for a C=C=O molecule, with AMS rdkit version
        img = classes[fmt.lower()]([rdmol], molsPerRow=1, subImgSize=size)
        img_text = img
        if isinstance(img, Image.Image):
            buf = BytesIO()
            img.save(buf, format=fmt)
            img_text = buf.getvalue()
    # If I do not make this correction to the SVG text, it is not readable in JupyterLab
    if fmt.lower() == "svg":
        img_text = _correct_svg(img_text)

    # Write to file, if required
    if filename is not None:
        mode = "w"
        if isinstance(img_text, bytes):
            mode = "wb"
        with open(filename, mode) as outfile:
            outfile.write(img_text)

    if as_string:
        img = img_text
    return img


@requires_optional_package("rdkit")
@requires_optional_package("PIL")
def get_reaction_image(
    reactants: Sequence[Molecule],
    products: Sequence[Molecule],
    filename: Optional[str] = None,
    fmt: str = "svg",
    size: Sequence[int] = (200, 100),
    as_string: bool = True,
):
    """
    Create a 2D reaction image from reactants and products (PLAMS molecules)

    * ``reactants`` -- Iterable of PLAMS Molecule objects representing the reactants.
    * ``products`` -- Iterable of PLAMS Molecule objects representing the products.
    * ``filename`` -- Optional: Name of image file to be created.
    * ``fmt``  -- The format of the image (svg, png, eps, pdf, jpeg).
                     The extension in the filename, if provided, takes precedence.
    * ``size`` -- Tuple/list containing width and height of image in pixels.
    * ``as_string`` -- Returns the image as a string or bytestring. If set to False, the original format
                       will be returned, which can be either a PIL image or SVG text
    *      Returns -- SVG image text file.
    """
    extensions = ["svg", "png", "eps", "pdf", "jpeg"]

    # Determine the type of image file
    if filename is not None:
        if "." in filename:
            extension = filename.split(".")[-1]
            if extension in extensions:
                fmt = extension
            else:
                msg = [f"Image type {extension} not available."]
                msg += [f"Available extensions are: {' '.join(extensions)}"]
                raise Exception("\n".join(msg))
        else:
            filename = ".".join([filename, fmt])

    if fmt.lower() not in extensions:
        raise Exception(f"Image type {fmt} not available.")

    # Get the actual image
    width = size[0]
    height = size[1]
    if fmt.lower() == "svg":
        img_text = _get_reaction_image_svg(reactants, products, width, height)
    else:
        img_text = _get_reaction_image_pil(reactants, products, fmt, width, height, as_string=as_string)

    # Write to file, if required
    if filename is not None:
        mode = "w"
        if isinstance(img_text, bytes):
            mode = "wb"
        with open(filename, mode) as outfile:
            outfile.write(img_text)
    return img_text


def _get_reaction_image_svg(
    reactants: Sequence[Molecule], products: Sequence[Molecule], width: int = 200, height: int = 100
):
    """
    Create a 2D reaction image from reactants and products (PLAMS molecules)

    * ``reactants`` -- Iterable of PLAMS Molecule objects representing the reactants.
    * ``products`` -- Iterable of PLAMS Molecule objects representing the products.
    *      Returns -- SVG image text file.
    """
    from rdkit import Chem

    def svg_arrow(x1, y1, x2, y2, prefix=""):
        """
        The reaction arrow in html format
        """
        # The arrow head
        l = ['<%sdefs> <%smarker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" ' % (prefix, prefix)]
        l += ['markerWidth="6" markerHeight="6" ']
        l += ['orient="auto-start-reverse"> ']
        l += ['<%spath d="M 0 0 L 10 5 L 0 10 z" /></%smarker></%sdefs>' % (prefix, prefix, prefix)]
        arrow = "".join(l)
        # The line
        l = ['<%sline x1="%i" y1="%i" x2="%i" y2="%i" ' % (prefix, x1, y1, x2, y2)]
        l += ['stroke="black" marker-end="url(#arrow)" />']
        line = "".join(l)
        return [arrow, line]

    def add_plus_signs_svg(img_text, width, height, nmols, nreactants, prefix=""):
        """
        Add the lines with + signs to the SVG image
        """
        y = int(0.55 * height)
        t = []
        for i in range(nmols - 1):
            x = int(((i + 1) * width) - (0.1 * width))
            if i + 1 in (nreactants, nreactants + 1):
                continue
            t += ['<%stext x="%i" y="%i" font-size="16">+</%stext>' % (prefix, x, y, prefix)]
        lines = img_text.split("\n")
        lines = lines[:-2] + t + lines[-2:]
        return "\n".join(lines)

    def add_arrow_svg(img_text, width, height, nreactants, prefix=""):
        """
        Add the arrow to the SVG image
        """
        y = int(0.5 * height)
        x1 = int((nreactants * width) + (0.3 * width))
        x2 = int((nreactants * width) + (0.7 * width))
        t = svg_arrow(x1, y, x2, y, prefix)
        lines = img_text.split("\n")
        lines = lines[:-2] + t + lines[-2:]
        return "\n".join(lines)

    # Get the rdkit molecules
    rdmols = [_rdmol_for_image(mol) for mol in reactants]
    rdmols += [Chem.Mol()]  # This is where the arrow will go
    rdmols += [_rdmol_for_image(mol) for mol in products]
    nmols = len(rdmols)

    # Place the molecules in a row of images
    subimg_size = [width, height]
    kwargs = {"legendFontSize": 16}  # ,"legendFraction":0.1}
    img_text = _MolsToGridSVG(rdmols, molsPerRow=nmols, subImgSize=subimg_size, **kwargs)
    img_text = _correct_svg(img_text)

    # Add + and =>
    nreactants = len(reactants)
    img_text = add_plus_signs_svg(img_text, width, height, nmols, nreactants)
    img_text = add_arrow_svg(img_text, width, height, nreactants)

    return img_text


def _get_reaction_image_pil(
    reactants: Sequence[Molecule],
    products: Sequence[Molecule],
    fmt: str,
    width: int = 200,
    height: int = 100,
    as_string: bool = True,
):
    """
    Create a 2D reaction image from reactants and products (PLAMS molecules)

    * ``reactants`` -- Iterable of PLAMS Molecule objects representing the reactants.
    * ``products`` -- Iterable of PLAMS Molecule objects representing the products.
    *      Returns -- SVG image text file.
    """
    from io import BytesIO
    from PIL import Image
    from PIL import ImageDraw
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage

    def add_arrow_pil(img, width, height, nreactants):
        """
        Add the arrow to the PIL image
        """
        y1 = int(0.5 * height)
        y2 = y1
        x1 = int((nreactants * width) + (0.3 * width))
        x2 = int((nreactants * width) + (0.7 * width))

        # Draw a line
        black = (0, 0, 0)
        draw = ImageDraw.Draw(img)
        draw.line(((x1, y1), (x2, y2)), fill=128)

        # Draw the arrow head
        headscale = 20
        xshift = width // headscale
        yshift = img.size[1] // headscale
        p1 = (x2, y2 + yshift)
        p2 = (x2, y2 - yshift)
        p3 = (x2 + xshift, y2)
        draw.polygon((p1, p2, p3), fill=black)

        return img

    def add_plus_signs_pil(img, width, height, nmols, nreactants):
        """
        Add the lines with + signs to the SVG image
        """
        black = (0, 0, 0)

        I1 = ImageDraw.Draw(img)
        y = int(0.5 * height)
        # myfont = ImageFont.truetype('FreeMono.ttf', 25)
        for i in range(nmols):
            x = int(((i + 1) * width) - (0.05 * width))
            if i + 1 in (nreactants, nreactants + 1):
                continue
            # I1.text((x, y), "+", font=myfont, fill=black)
            I1.text((x, y), "+", fill=black)
        return img

    def join_pil_images(pil_images):
        """
        Create a new image which connects the ones above with text
        """
        white = (255, 255, 255)

        widths = [img.width for img in pil_images]
        width = sum(widths)
        height = max([img.height for img in pil_images])
        final_img = Image.new("RGB", (width, height), white)

        # Concatenate the PIL images
        for i, img in enumerate(pil_images):
            pos = sum(widths[:i])
            h = int((height - img.height) / 2)
            final_img.paste(img, (pos, h))

        return final_img

    nreactants = len(reactants)
    nmols = nreactants + len(products)

    if not hasattr(rdMolDraw2D, "MolDraw2DCairo"):
        # We are working with the old AMS version of RDKit
        white = (255, 255, 255)
        rimages = [to_image(mol, fmt=fmt, as_string=False) for i, mol in enumerate(reactants)]
        pimages = [to_image(mol, fmt=fmt, as_string=False) for i, mol in enumerate(products)]
        blanc = Image.new("RGB", (width, height), white)

        # Get the image (with arrow)
        nreactants = len(reactants)
        all_images = rimages + [blanc] + pimages
        img = join_pil_images(all_images)

    else:
        # We have a later version of RDKit that can regulate the font sizes
        rdmols = [_rdmol_for_image(mol) for mol in reactants]
        rdmols += [Chem.Mol()]  # This is where the arrow will go
        rdmols += [_rdmol_for_image(mol) for mol in products]

        # Place the molecules in a row of images
        subimg_size = [width, height]
        kwargs = {"legendFontSize": 16}  # ,"legendFraction":0.1}
        img = MolsToGridImage(rdmols, molsPerRow=nmols + 1, subImgSize=subimg_size, **kwargs)

    # Add + and =>
    img = add_plus_signs_pil(img, width, height, nmols, nreactants)
    img = add_arrow_pil(img, width, height, nreactants)

    # Get the bytestring
    img_text = img
    if as_string:
        buf = BytesIO()
        img.save(buf, format=fmt)
        img_text = buf.getvalue()

    return img_text


def _correct_svg(image):
    """
    Correct for a bug in the AMS rdkit created SVG file
    """
    if not "svg:" in image:
        return image
    image = image.replace("svg:", "")
    lines = image.split("\n")
    for iline, line in enumerate(lines):
        if "xmlns:svg=" in line:
            lines[iline] = line.replace("xmlns:svg", "xmlns")
            break
    image = "\n".join(lines)
    return image


def _presanitize(mol, rdmol):
    """
    Change bonding and atom charges to avoid failed sanitization

    Note: Used by to_rdmol
    """
    from rdkit import Chem

    mol = mol.copy()
    for i in range(10):
        try:
            rdmol_test = copy.deepcopy(rdmol)
            Chem.SanitizeMol(rdmol_test)
            stored_exc = None
            break
        except ValueError as exc:
            stored_exc = exc
            text = repr(exc)
        # Fix the problem
        rdmol, bonds, charges = _kekulize(mol, rdmol, text)
        # print ("REB bonds charges: ",bonds, charges)
        # print (rdmol.Debug())
        # rdmol = to_rdmol(mol, sanitize=False)
        # print (rdmol.Debug())

    if stored_exc is not None:
        raise stored_exc
    else:
        rdmol = rdmol_test
    return rdmol


def _kekulize(mol, rdmol, text, use_dfs=True):
    """
    Kekulize the atoms indicated as problematic by RDKit

    * ``mol`` - PLAMS molecule
    * ``text`` - Sanitation error produced by RDKit

    Note: Returns the changes in bond orders and atomic charges, that will make sanitation succeed.
    """
    from rdkit import Chem

    # Find the indices of the problematic atoms
    indices = _find_aromatic_sequence(rdmol, text)
    if indices is None:
        return rdmol, {}, {}

    # Set the bond orders along the chain to 2, 1, 2, 1,...
    altered_bonds = {}
    if len(indices) > 1:
        emol = Chem.RWMol(rdmol)
        if use_dfs:
            # It may be better to create the tree first, and then start at a leaf
            altered_bonds = _alter_aromatic_bonds(emol, indices[0])
            indices = list(set([ind for tup in altered_bonds.keys() for ind in tup]))
        else:
            indices = _order_atom_indices(emol, indices)
            altered_bonds = _alter_bonds_along_chain(emol, indices)
        _adjust_atom_aromaticity(emol, altered_bonds)
        rdmol = emol.GetMol()
        # Adjust the PLAMS molecule as well
        for (iat, jat), order in altered_bonds.items():
            for bond in mol.atoms[iat].bonds:
                if mol.index(bond.other_end(mol.atoms[iat])) - 1 == jat:
                    bond.order = order
                    break

    # If the atom at the chain end has the wrong bond order, give it a charge.
    altered_charge = {}
    atom_indices, dangling_bonds = _get_charged_atoms(rdmol, indices)
    if len(atom_indices) > 0:
        iat = atom_indices[0]
        ndangling = dangling_bonds[iat]
        altered_charge = _guess_atomic_charge(iat, ndangling, rdmol, mol)
        rdmol.GetAtomWithIdx(iat).SetFormalCharge(altered_charge[iat])
        for k, v in altered_charge.items():
            mol.atoms[k].properties.rdkit.charge = v

    return rdmol, altered_bonds, altered_charge


def _find_aromatic_sequence(rdmol, text):
    """
    Find the sequence of atoms with 1.5 bond orders
    """
    indices = None
    lines = text.split("\n")
    line = lines[-1]
    if "Unkekulized atoms:" in line:
        text = line.split("Unkekulized atoms:")[-1].split("\\n")[0]
        if '"' in text:
            text = text.split('"')[0]
        indices = [int(w) for w in text.split()]
        if len(indices) > 1:
            return indices
        line = "atom %i marked aromatic" % (indices[0])
    iat: Any
    if "marked aromatic" in line:
        iat = int(line.split("atom")[-1].split()[0])
        indices = [iat]
        while iat is not None:
            icurrent = iat
            iat = None
            for bond in rdmol.GetAtomWithIdx(icurrent).GetBonds():
                if str(bond.GetBondType()) == "AROMATIC":
                    iat = bond.GetOtherAtomIdx(icurrent)
                    if iat in indices:
                        iat = None
                        continue
                    indices.append(iat)
                    break
    elif "Explicit valence for atom" in line:
        iat = int(line.split("#")[-1].split()[0])
        indices = [iat]
    return indices


def _alter_aromatic_bonds(emol, iat, depth=0, double_first=False):
    """
    Switch all thearomitic bonds to single/double, starting at iat

    * ``emol`` -- RDKit EditableMol type, for which bond orders will be changed
    * ``iat``  -- Starting point for depth first search
    """
    from collections import OrderedDict
    from rdkit import Chem
    from scm.plams import PeriodicTable as PT

    # Use OrderedDict, so that the leaves of the tree
    # will be at the end
    bonds_changed = OrderedDict()

    at = emol.GetAtomWithIdx(iat)
    valence = PT.get_connectors(at.GetAtomicNum())
    bonds = at.GetBonds()
    are_aromatic = [str(b.GetBondType()) == "AROMATIC" for b in bonds]
    int_orders = sum([b.GetBondTypeAsDouble() for i, b in enumerate(bonds) if not are_aromatic[i]])
    numbonds = len([b for i, b in enumerate(bonds) if are_aromatic[i]])
    valence -= int(int_orders)
    # Here I place the double bond first.
    # I could also start with a single bond instead.
    orders = [1 for i in range(numbonds)]
    if valence > numbonds:
        for i in range(min(numbonds, valence - numbonds)):
            if double_first:
                orders[i] = 2
            else:
                orders[numbonds - i - 1] = 2

    for i, bond in enumerate(bonds):
        jat = bond.GetOtherAtomIdx(iat)
        if are_aromatic[i]:
            pair = tuple(sorted([iat, jat]))
            order = orders.pop(0)
            bond.SetBondType(Chem.BondType(order))
            bond.SetIsAromatic(False)  # This is necessary with the newer RDKit
            bonds_changed[pair] = order
            d = _alter_aromatic_bonds(emol, jat, depth + 1)
            bonds_changed.update(d)
    return bonds_changed


def _alter_bonds_along_chain(emol, indices):
    """
    Along the chain of atoms (indices), alternate double and single bonds
    """
    from scm.plams import PeriodicTable as PT
    from rdkit import Chem

    if len(indices) > 1:
        # The first bond order is set to 2, if aromic
        # Else it is flipped
        first_bond = emol.GetBondBetweenAtoms(indices[0], indices[1])
        if str(first_bond.GetBondType()) == "AROMATIC":
            new_order = 2
        else:
            new_order = first_bond.GetBondTypeAsDouble()
            new_order = (new_order % 2) + 1
        # Unless this breaks valence rules.
        iat = indices[0]
        at = emol.GetAtomWithIdx(iat)
        valence = PT.get_connectors(at.GetAtomicNum())
        orders = [b.GetBondTypeAsDouble() for b in at.GetBonds()]
        jat = indices[1]
        at_next = emol.GetAtomWithIdx(jat)
        valence_next = PT.get_connectors(at_next.GetAtomicNum())
        orders_next = [b.GetBondTypeAsDouble() for b in at_next.GetBonds()]
        if sum(orders) > valence or sum(orders_next) > valence_next:
            new_order = 1
    # Set the bond orders along the chain to 2, 1, 2, 1,...
    altered_bonds = {}
    for i, iat in enumerate(indices[:-1]):
        bond = emol.GetBondBetweenAtoms(iat, indices[i + 1])
        bond.SetBondType(Chem.BondType(new_order))
        bond.SetIsAromatic(False)  # This is necessary with the newer RDKit versions
        altered_bonds[iat, indices[i + 1]] = new_order
        new_order = ((new_order) % 2) + 1
    return altered_bonds


def _order_atom_indices(rdmol, indices):
    """
    Order the atomic indices so that they are consecutive along a bonded chain
    """
    # Order the indices, so that they are consecutive in the molecule
    start = 0
    for i, iat in enumerate(indices):
        at = rdmol.GetAtomWithIdx(iat)
        neighbors = [b.GetOtherAtomIdx(iat) for b in at.GetBonds()]
        relevant_neighbors = [jat for jat in neighbors if jat in indices]
        if len(relevant_neighbors) == 1:
            start = i
            break
    iat = indices[start]
    atoms = [iat]
    while 1:
        at = rdmol.GetAtomWithIdx(iat)
        neighbors = [b.GetOtherAtomIdx(iat) for b in at.GetBonds()]
        relevant_neighbors = [jat for jat in neighbors if jat in indices and not jat in atoms]
        if len(relevant_neighbors) == 0:
            break
        iat = relevant_neighbors[0]
        atoms.append(iat)
    if len(atoms) < len(indices):
        raise Exception("The unkekulized atoms are not in a consecutive chain")
    indices = atoms
    return indices


def _adjust_atom_aromaticity(emol, altered_bonds):
    """
    Assign aromaticity to the atoms based on the new bond orders
    """
    indices = list(set([ind for tup in altered_bonds.keys() for ind in tup]))
    for ind in indices:
        at = emol.GetAtomWithIdx(ind)
        if not at.GetIsAromatic():
            continue
        # Only change this if we are sure the atom is not aromatic anymore
        aromatic = False
        for bond in at.GetBonds():
            if str(bond.GetBondType()) == "AROMATIC":
                aromatic = True
                break
        if not aromatic:
            at.SetIsAromatic(False)


def _get_charged_atoms(rdmol, indices):
    """
    Locate the atoms that need a charge
    """
    from scm.plams import PeriodicTable as PT

    atom_indices = []
    dangling_bonds = {}
    for iat in indices[::-1]:
        at = rdmol.GetAtomWithIdx(iat)
        valence = PT.get_connectors(at.GetAtomicNum())
        bonds = at.GetBonds()
        orders = [b.GetBondTypeAsDouble() for b in bonds]
        ndangling = int(valence - sum(orders))
        if ndangling != 0:
            dangling_bonds[iat] = ndangling
            atom_indices.append(iat)
    return atom_indices, dangling_bonds


def _guess_atomic_charge(iat, ndangling, rdmol, mol):
    """
    Guess the best atomic charge for atom iat
    """
    # Get the total charge from atomic charges already set
    charges = [at.GetFormalCharge() for at in rdmol.GetAtoms()]
    totcharge = sum(charges) - charges[iat]

    # Here I hope that the estimated charge will be more reliable than the
    # actual (user defined) system charge, but am not sure
    est_charges = mol.guess_atomic_charges(adjust_to_systemcharge=False, depth=0)
    molcharge = int(sum(est_charges))
    molcharge = molcharge - totcharge

    # Adjust the sign based on the estimated charge
    sign = ndangling / (abs(ndangling))
    if molcharge != 0:
        sign = molcharge / (abs(molcharge))
    elif est_charges[iat] != 0:
        sign = est_charges[iat] / abs(est_charges[iat])

    # Then use the over/under valence with the new sign as the atomic charge
    charge = int(sign * abs(ndangling))
    # Perhaps we tried this already
    if charge == charges[iat]:
        charge = -charge
    altered_charge = {iat: charge}
    return altered_charge


def _rdmol_for_image(mol, remove_hydrogens=True):
    """
    Convert PLAMS molecule to an RDKit molecule specifically for a 2D image
    """
    from rdkit.Chem import AllChem
    from rdkit.Chem import RemoveHs

    rdmol = to_rdmol(mol, presanitize=True)

    # Flatten the molecule
    AllChem.Compute2DCoords(rdmol)
    # Remove the Hs only if there are carbon atoms in this system
    # Otherwise this will turn an OH radical into a water molecule.
    carbons = [i for i, at in enumerate(mol.atoms) if at.symbol in ["C", "Si"]]
    if remove_hydrogens and len(carbons) > 0:
        rdmol = RemoveHs(rdmol)
    else:
        for atom in rdmol.GetAtoms():
            atom.SetNoImplicit(True)

    ids = [c.GetId() for c in rdmol.GetConformers()]
    for cid in ids:
        rdmol.RemoveConformer(cid)
    return rdmol


def _MolsToGridSVG(
    mols,
    molsPerRow=3,
    subImgSize=(200, 200),
    legends=None,
    highlightAtomLists=None,
    highlightBondLists=None,
    drawOptions=None,
    **kwargs,
):
    """
    Replaces the old version of this function in our RDKit for a more recent one, with more options
    """
    from rdkit.Chem.Draw import rdMolDraw2D

    if legends is None:
        legends = [""] * len(mols)

    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow:
        nRows += 1

    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

    d2d = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
    if drawOptions is not None:
        d2d.SetDrawOptions(drawOptions)
    else:
        dops = d2d.drawOptions()
        for k, v in list(kwargs.items()):
            if hasattr(dops, k):
                setattr(dops, k, v)
                del kwargs[k]

    d2d.DrawMolecules(
        list(mols),
        legends=legends or None,
        highlightAtoms=highlightAtomLists or [],
        highlightBonds=highlightBondLists or [],
        **kwargs,
    )
    d2d.FinishDrawing()
    res = d2d.GetDrawingText()
    return res
