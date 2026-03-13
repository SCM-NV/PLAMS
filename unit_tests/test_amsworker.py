from typing import Any, Dict
from unittest.mock import MagicMock
from abc import ABC, abstractmethod
import numpy as np
import pytest

from scm.plams.interfaces.adfsuite.amsworker import AMSWorker, AMSWorkerResults
from scm.plams.mol.molecule import Molecule
from scm.plams.tools.units import Units
from scm.plams.interfaces.molecule.ase import fromASE
from test_helpers import skip_if_no_scm_base


class PrepareSystemTestBase(ABC):

    @abstractmethod
    @pytest.fixture
    def mol(self, folder): ...

    @property
    def expected_atoms_symbols(self):
        return []

    @property
    def expected_coords(self):
        return []

    @property
    def expected_charge(self):
        return 0.0

    @property
    def expected_atomic_info(self):
        return None

    @property
    def expected_lattice_vectors(self):
        return None

    @property
    def expected_bonds(self):
        return None

    @property
    def expected_bond_orders(self):
        return None

    def _prepare_and_capture_system_payload(self, system) -> Dict[str, Any]:
        worker = AMSWorker.__new__(AMSWorker)
        worker._prune_restart_cache = MagicMock()
        worker._call = MagicMock()

        worker._prepare_system(system)

        worker._prune_restart_cache.assert_called_once_with()
        worker._call.assert_called_once()
        request, payload = worker._call.call_args.args

        assert request == "SetSystem"
        assert isinstance(payload, dict)
        return payload

    def test_prepare_system_with_molecule(self, mol) -> None:
        payload = self._prepare_and_capture_system_payload(mol)

        assert np.all(payload["atomSymbols"] == self.expected_atoms_symbols)
        assert np.allclose(payload["coords"], self.expected_coords)
        assert payload["totalCharge"] == self.expected_charge
        if self.expected_atomic_info is not None:
            assert np.all(payload["atomicInfo"] == self.expected_atomic_info)
        else:
            assert "atomicInfo" not in payload
        if self.expected_lattice_vectors is not None:
            assert np.allclose(payload["latticeVectors"], self.expected_lattice_vectors)
        else:
            assert "latticeVectors" not in payload

        if self.expected_bonds is not None:
            assert np.all(payload["bonds"] == self.expected_bonds)
        else:
            assert "bonds" not in payload
        if self.expected_bond_orders is not None:
            assert np.all(payload["bondOrders"] == self.expected_bond_orders)
        else:
            assert "bondOrders" not in payload


class TestAMSWorkerWithHydroxideMolecule(PrepareSystemTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "hydroxide.xyz")
        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]]
        mol.properties.charge = -1
        mol[1].properties.region = {"oxygen", "hydroxide"}
        mol[2].properties.region = {"hydrogen", "hydroxide"}
        mol[2].properties.mass = 2.0141
        mol[1].properties.adf.f = "myfrag"
        mol[2].properties.adf.f = "myfrag"
        mol.guess_bonds()
        return mol

    @property
    def expected_atoms_symbols(self):
        return ["O", "H"]

    @property
    def expected_coords(self):
        return np.array([[1.88972612, 0, 0], [0, 0, 0]])

    @property
    def expected_charge(self):
        return -1

    @property
    def expected_atomic_info(self):
        return ["region=hydroxide,oxygen adf.f=myfrag", "region=hydrogen,hydroxide mass=2.0141 adf.f=myfrag"]

    @property
    def expected_lattice_vectors(self):
        return np.array([[9.44863062, 0, 0], [0, 11.33835675, 0], [0, 0, 13.22808287]])

    @property
    def expected_bonds(self):
        return [[1, 2]]

    @property
    def expected_bond_orders(self):
        return [[1]]


class TestAMSWorkerWithHydroxideChemicalSystem(TestAMSWorkerWithHydroxideMolecule):

    @pytest.fixture
    def mol(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "hydroxide.xyz"))

        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]]
        mol.charge = -1
        mol.add_atoms_to_region([0, 1], "hydroxide")
        mol.add_atom_to_region(0, "oxygen")
        mol.add_atom_to_region(1, "hydrogen")
        mol.enable_atom_attributes("adf")
        mol.atoms[1].mass = 2.0141
        mol.atoms[0].adf.f = "myfrag"
        mol.atoms[1].adf.f = "myfrag"
        mol.guess_bonds()
        return mol

    @property
    def expected_atomic_info(self):
        return ["adf.f=myfrag region=hydroxide,oxygen", "mass=2.0141 adf.f=myfrag region=hydrogen,hydroxide"]


class TestAMSWorkerWithWaterMolecule(PrepareSystemTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "water.xyz")
        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0]]
        mol[1].properties.region = {"oxygen"}
        mol[2].properties.mass = 2
        mol[1].properties.forcefield.type = "o"
        mol[2].properties.forcefield.type = "h2"
        mol.guess_bonds()
        return mol

    @property
    def expected_atoms_symbols(self):
        return ["O", "H", "H"]

    @property
    def expected_coords(self):
        return np.array([[0.0, 0.0, 0.0], [1.88972612, 0.0, 0.0], [0.0, 1.88972612, 0.0]])

    @property
    def expected_charge(self):
        return 0

    @property
    def expected_atomic_info(self):
        return ["region=oxygen forcefield.type=o", "mass=2 forcefield.type=h2", ""]

    @property
    def expected_lattice_vectors(self):
        return np.array([[9.44863062, 0.0, 0.0], [0.0, 11.33835675, 0.0]])

    @property
    def expected_bonds(self):
        return [[1, 3], [1, 2]]

    @property
    def expected_bond_orders(self):
        return [1, 1]


class TestAMSWorkerWithWaterChemicalSystem(TestAMSWorkerWithWaterMolecule):

    @pytest.fixture
    def mol(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "water.xyz"))

        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0]]
        mol.add_atom_to_region(0, "oxygen")
        mol.enable_atom_attributes("adf")
        mol.enable_atom_attributes("forcefield")
        mol.atoms[1].mass = 2
        mol.atoms[0].forcefield.type = "o"
        mol.atoms[1].forcefield.type = "h2"
        mol.guess_bonds()
        return mol

    @property
    def expected_atomic_info(self):
        return ["forcefield.type=o region=oxygen", "mass=2 forcefield.type=h2", ""]

    @property
    def expected_bonds(self):
        return [[1, 2], [1, 3]]


class TestAMSWorkerWithBenzeneMolecule(PrepareSystemTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "benzene.xyz")
        mol.guess_bonds()
        return mol

    @property
    def expected_atoms_symbols(self):
        return ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]

    @property
    def expected_coords(self):
        return np.array(
            [
                [
                    [2.25606881, -1.30254194, 0.0],
                    [2.25606881, 1.30254194, 0.0],
                    [0.0, 2.60508388, 0.0],
                    [-2.25606881, 1.30254194, 0.0],
                    [-2.25606881, -1.30254194, 0.0],
                    [-0.0, -2.60508388, 0.0],
                    [4.03061813, -2.32707846, -0.0],
                    [4.03061813, 2.32707846, -0.0],
                    [0.0, 4.65415692, -0.0],
                    [-4.03061813, 2.32707846, -0.0],
                    [-4.03061813, -2.32707846, -0.0],
                    [-0.0, -4.65415692, -0.0],
                ]
            ]
        )

    @property
    def expected_charge(self):
        return 0

    @property
    def expected_bonds(self):
        return np.array(
            [[3, 4], [5, 6], [1, 6], [2, 3], [4, 5], [1, 2], [3, 9], [6, 12], [5, 11], [4, 10], [2, 8], [1, 7]]
        )

    @property
    def expected_bond_orders(self):
        return np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])


class TestAMSWorkerWithBenzeneChemicalSystem(TestAMSWorkerWithBenzeneMolecule):

    @pytest.fixture
    def mol(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "benzene.xyz"))
        mol.guess_bonds()
        return mol

    @property
    def expected_bonds(self):
        return np.array(
            [[1, 2], [1, 6], [1, 7], [2, 3], [2, 8], [3, 4], [3, 9], [4, 5], [4, 10], [5, 6], [5, 11], [6, 12]]
        )

    @property
    def expected_bond_orders(self):
        return np.array([1.5, 1.5, 1.0, 1.5, 1.0, 1.5, 1.0, 1.5, 1.0, 1.5, 1.0, 1.0])


class TestAMSWorkerWithChlorophyllMolecule(PrepareSystemTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "chlorophyl1.xyz")
        return mol

    @property
    def expected_atoms_symbols(self):
        return (
            ["H", "H", "H", "H", "O", "H", "C", "C", "C", "C", "H", "H", "H"]
            + ["C", "H", "H", "H", "C", "C", "H", "C", "C", "C", "C", "C", "H"]
            + ["C", "N", "C", "H", "N", "H", "H", "H", "C", "C", "H", "H", "H"]
            + ["H", "C", "C", "H", "H", "C", "H", "H", "H", "C", "C", "H", "C"]
            + ["O", "C", "C", "H", "O", "H", "H", "C", "C", "C", "H", "H", "H"]
            + ["H", "C", "H", "C", "C", "H", "C", "H", "N", "H", "H", "C", "H"]
            + ["C", "H", "H", "C", "C", "C", "H", "H", "H", "H", "H", "O", "C"]
            + ["H", "H", "H", "C", "C", "N", "H", "H", "H", "Mg", "H", "H", "C"]
            + ["H", "C", "C", "O", "H", "H", "H", "H", "C", "C", "H", "C", "C"]
            + ["C", "H", "H", "H", "C", "H", "H", "C", "C", "H", "H", "C", "C"]
            + ["C", "H", "H", "C", "C", "H", "H"]
        )

    @property
    def expected_coords(self):
        return np.array(
            [
                [
                    [-20.346418511914344, 15.807621393456682, -11.424590619999668],
                    [-0.5418676278796919, -5.518081542130608, 8.538205929711145],
                    [-14.012882352485224, 1.150217710549843, 2.5212839338324504],
                    [-20.473564954757542, 16.93172119923807, -4.6734703404552445],
                    [0.1098478898983714, -11.688403945901925, 7.6873283962065235],
                    [-22.981097351581088, 15.99125930906944, -9.431881980486294],
                    [4.316030533708405, 1.3804241470517544, -4.233283206163291],
                    [6.423009022251776, -3.3117223566931666, 9.764488896229425],
                    [6.430363836328819, -4.339079523250465, 6.936734848683541],
                    [2.881849347589421, 0.23038218103150152, -0.6865715161467855],
                    [18.532277652821353, -2.7730124611677245, 6.504340944929545],
                    [2.459399182703206, -7.154743103050993, 8.319332180778616],
                    [-7.929135861387512, -11.367192078690286, 1.310780180454796],
                    [2.2495564349202617, -0.9632708844928647, 1.8284234091429101],
                    [-8.484851402308461, -8.518913715704842, 7.918836854008302],
                    [7.982006618902292, -2.3725908337162718, 10.453524615242724],
                    [-12.002406507948116, -2.0393017268421487, -1.9174597452307778],
                    [8.824971869123106, 1.6968039845627252, -6.045400168576806],
                    [-13.881886537526166, 3.902195630224358, -4.220684402090412],
                    [16.026795608843024, 5.9493189334963335, -7.5463888423796845],
                    [-16.0497690093401, 6.083610430816739, -4.1473082263973176],
                    [-0.23113807148135185, 1.695345115994514, -7.189623667859336],
                    [-4.7801642774516155, -10.521609557786865, 7.120559847182638],
                    [-9.864213523278176, -4.090909350207862, 3.3070490639869674],
                    [-14.37727824109681, 1.8114763452572662, -2.268121104368585],
                    [-4.671310383494797, -10.050200368463845, 2.1018535512933867],
                    [18.406782830611082, -0.8151824967365662, 1.5976991876828515],
                    [5.160832597722355, 0.6193501784415977, -1.6079018190297059],
                    [6.245372876810829, 1.6793089001009398, -6.25663186506135],
                    [-9.374292576838323, -3.71406016643499, 1.3381320763826294],
                    [11.76461093132971, -1.793737486125404, 1.171900428103799],
                    [-3.3520774580845503, 6.94106043178281, 4.351375971143406],
                    [-15.553611966688983, 2.5829740431708585, -0.922443351570325],
                    [-16.148416492497816, 12.077831146760305, -8.366615018142877],
                    [1.611583005520477, 1.0657601808619435, -4.786766980531058],
                    [-18.861088108780116, 4.88726828501703, -3.0197483320817375],
                    [-14.290823270895181, 8.747066019909283, -6.434054471450215],
                    [9.053332043201257, 6.644872317913465, -9.934176853590197],
                    [-19.931551817966252, 6.592708317969295, -3.0537369461592565],
                    [-1.9127089737534688, 2.138204323026688, -6.193216435771157],
                    [12.517971037899143, 0.778608737320559, -4.369403958372335],
                    [12.893680716818864, 2.8560375784531575, -6.00791178171648],
                    [-1.5913572668086073, -2.9393234779437893, 2.0319865970137228],
                    [12.187093221560042, 7.156358818887549, -8.234848194924972],
                    [14.977381558959095, 4.967916126837684, -5.972424614822132],
                    [-13.555651778275271, -4.144959296824407, 3.8309189402557235],
                    [13.82974899374774, 6.589202876008114, -5.023256756397349],
                    [-18.94779441255632, 10.977456852473592, -10.193428380627605],
                    [-12.485371372523224, -2.4924051391959194, 3.6326942286869786],
                    [-1.8658003021618832, 1.0425996974785299, 3.777065525156138],
                    [5.040241614805486, 2.2709481346510256, -7.905920575946485],
                    [10.826218291267542, 5.783276257829965, -8.98699698704852],
                    [-3.315911879511462, -1.7819947279869797, -2.1125172758146498],
                    [19.064016348673054, 1.50603046310338, 0.3551343408747964],
                    [12.142106401437202, 4.2121220530197325, -11.280591599577061],
                    [-17.62022590521035, 15.617675572039921, -7.787247665046111],
                    [-2.4749365108998784, 3.3033206343430805, 2.7795868183552974],
                    [-17.121925913133904, -2.519246809070104, -1.0542101747882306],
                    [-20.666925511281498, 18.293467293274023, -9.143747879360104],
                    [-0.616748025567988, -9.626374482958902, 6.9560081654285995],
                    [-19.282903325688455, 14.417555863717086, -7.531433549829396],
                    [-0.794404957996306, -0.9261018613476004, 1.594050126536199],
                    [-20.306976148241155, 9.698520446944864, -5.347049989495228],
                    [4.583502369387936, -2.7436631247261616, 10.485786019642335],
                    [0.10505932389856969, -5.705081280519075, 3.9280810986794816],
                    [-13.329615967330408, 2.8208735547261212, -6.025686545644709],
                    [-20.699838871194103, 14.895325871175595, -5.116970164643667],
                    [2.6520038487773134, -7.768180328448637, 4.071157932753273],
                    [14.113296729569464, -0.8836170386137638, 0.07766207454374528],
                    [-6.9797922584331396, -8.899992106544625, 6.338417322009029],
                    [-0.4653186020233512, 0.16833302372941436, -8.384087425891167],
                    [4.142105810376223, -1.7027018197976824, 3.569601942564098],
                    [13.422068928251601, 2.7393564388881395, -11.173325075564929],
                    [6.654814167055121, -1.5973609267065434, 2.9393064704086678],
                    [-19.353316410818138, 3.173634399588388, -4.326767957548528],
                    [-12.170482528924582, 4.765016232836615, -3.8617082480042506],
                    [16.753792146247804, -3.9223816846876107, 6.323977924690763],
                    [7.183761737420747, -6.440315134101099, 6.471060978148761],
                    [14.365597843920499, 0.18952630221709235, -2.2450588867436525],
                    [-11.832605167019867, -1.1127765668426326, 5.191321439017067],
                    [-10.216105094123115, 0.7624024460757686, -0.810174722506308],
                    [-1.4788788781447568, -0.5305443889409343, -1.3268900956672307],
                    [-10.038257299356209, -6.902736785975399, 3.8512278269170763],
                    [-18.316317861572998, 9.919087390485059, -5.743652480448812],
                    [-6.637540180549917, -7.91144613513489, 0.09402710278300445],
                    [-22.90494894766317, 14.978386893257397, -5.598621669562158],
                    [-16.47537879630669, 6.889380201726545, -8.051489958778658],
                    [-18.967118751906742, 4.467936168236447, -0.839023281524845],
                    [-11.38946628052125, -7.5054422567112935, 2.5521809559700817],
                    [-2.3893980578686933, 0.14253448267602337, 5.7361540471862975],
                    [-7.690943442282933, -8.482557274793166, 4.006325208869612],
                    [-20.006534260865276, 13.890739904420037, -3.527262628741857],
                    [20.10107158970193, -1.9220820153316127, 2.1955782978964513],
                    [-8.291613677982605, -3.1317958631891734, 4.088751282973538],
                    [-12.33823729645986, 0.015578902171414848, -1.1424169211073878],
                    [10.73855389634615, 3.7155265939077, -7.0639908956419655],
                    [10.029619361240545, 0.08711826407137263, -4.333011612664292],
                    [-0.08998497860242992, 5.924550293180863, 5.323772343092088],
                    [-11.314217386238651, -6.885906885109482, 5.460484579578139],
                    [-15.967781351697088, 0.8673502761329852, -2.861060470492413],
                    [8.68341248512739, -0.5095306344350956, -0.424785866376253],
                    [-16.855003987756636, -4.082547312106392, 2.1404946710897343],
                    [15.323927094597469, -3.831463181379615, 7.642180949363089],
                    [-18.285808233290915, 11.82416187069504, -8.338854941372125],
                    [10.806077590231281, 2.823008945246949, -12.179708171865004],
                    [0.6204386606893821, -6.813290379318487, 7.408990635310394],
                    [3.6103746734290234, -3.619728817745921, 5.7301862920847295],
                    [-2.88145250510325, -9.236695948525947, 5.74265092560276],
                    [20.81884248416629, 2.232478979932018, 0.1373396255594271],
                    [17.03038002075253, -6.001227820413678, 5.80668240560958],
                    [-4.359214555108353, -10.555324161576314, 9.243500852300723],
                    [12.667229166126585, 6.03338717960256, -12.42417826115559],
                    [10.478877180930702, -2.96169027635486, 5.320697758687322],
                    [1.9696747677803128, -5.971833130545124, 4.884496056792204],
                    [-15.226320850534425, 7.26218726131884, -2.6732179142523624],
                    [-13.961143540919096, -0.9047347280564568, 1.4880516087536904],
                    [-6.650237250381277, -9.457346040015624, 1.6115943438772247],
                    [-20.973095158541113, 16.21077745461047, -9.499760942882851],
                    [11.265008917953395, -4.238067992710843, 7.070300691172091],
                    [-18.294398928253464, -0.9888086433410574, 1.3317183459156496],
                    [-2.2112290096606015, 4.672034148600529, 6.712456483034581],
                    [8.008047044899635, -2.8538379372440934, 4.964798078732051],
                    [-17.429647025541836, 11.179774710828276, -4.274515140476502],
                    [16.614114929472215, 0.22935228029358048, -2.323688501063206],
                    [0.7615275028800667, 0.4356215559748572, -2.3066469508713308],
                    [-2.1277881526277507, 5.085390951375045, 4.522796807360458],
                    [0.36965688614254544, 3.227058846857683, -8.286456615388502],
                    [6.131909940836049, -4.795907585795873, 11.15904441524113],
                    [15.084687767219846, -2.6373376843240925, 3.9566197426135803],
                    [-16.738819846162396, -2.0875350964470973, 0.8904994211596509],
                    [-16.428305718542266, 7.709037569926224, -6.319894226535447],
                    [2.407420375919249, -2.9600726707921803, 6.997396947010153],
                    [-4.953506965137413, -12.744023856379126, 6.351826818589372],
                    [15.84764012359788, -1.4444858626765678, 1.9094416372839906],
                    [12.352058753335374, -2.520972129021887, 3.538960033453291],
                    [17.65046908210886, 2.3813799499219064, -0.9387611366564685],
                    [16.827138086322908, 4.538511769812845, -5.340341461752806],
                ]
            ]
        )


class TestAMSWorkerWithChlorophyllChemicalSystem(TestAMSWorkerWithChlorophyllMolecule):

    @pytest.fixture
    def mol(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "chlorophyl1.xyz"))
        return mol


class WorkerResultsTestBase(ABC):

    @property
    def coord_shift(self):
        return np.array([0.1, -0.2, 0.3])

    @property
    def lattice_scaling(self):
        return 1.042

    @abstractmethod
    @pytest.fixture
    def mol(self, xyz_folder): ...

    @abstractmethod
    @pytest.fixture
    def system(self, xyz_folder): ...

    def _create_results_from_molecule(self, molecule, include_xyz=False) -> AMSWorkerResults:
        results = {}

        if include_xyz:
            shifted_coords = molecule.as_array() + self.coord_shift
            results["xyzAtoms"] = shifted_coords * Units.conversion_ratio("Angstrom", "au")
            if molecule.lattice:
                scaled_lattice = np.asarray(molecule.lattice) * self.lattice_scaling
                results["latticeVectors"] = scaled_lattice * Units.conversion_ratio("Angstrom", "Bohr")

        return AMSWorkerResults("test", molecule, results)

    def _create_results_from_chemical_system(self, molecule, include_xyz=False) -> AMSWorkerResults:
        skip_if_no_scm_base()

        results = {}

        if include_xyz:
            shifted_coords = molecule.coords + self.coord_shift
            results["xyzAtoms"] = shifted_coords * Units.conversion_ratio("Angstrom", "au")
            if molecule.has_lattice():
                scaled_lattice = np.asarray(molecule.lattice.vectors) * self.lattice_scaling
                results["latticeVectors"] = scaled_lattice * Units.conversion_ratio("Angstrom", "Bohr")
            else:
                results["latticeVectors"] = np.zeros((0, 3))

        return AMSWorkerResults("test", molecule, results)

    @pytest.mark.parametrize("include_xyz", [True, False])
    def test_get_input_molecule(self, mol, system, include_xyz):
        results_mol = self._create_results_from_molecule(mol, include_xyz=include_xyz)
        results_sys = self._create_results_from_chemical_system(system, include_xyz=include_xyz)

        input_mol_from_mol = results_mol.get_input_molecule()
        input_mol_from_sys = results_sys.get_input_molecule()
        for input_mol in [input_mol_from_mol, input_mol_from_sys]:
            assert np.all(input_mol.symbols == mol.symbols)
            assert np.allclose(input_mol.as_array(), mol.as_array())
            assert np.allclose(input_mol.lattice, mol.lattice)
            assert input_mol.label(4) == mol.label(4)

    @pytest.mark.parametrize("include_xyz", [True, False])
    def test_get_input_system(self, mol, system, include_xyz):
        results_mol = self._create_results_from_molecule(mol, include_xyz=include_xyz)
        results_sys = self._create_results_from_chemical_system(system, include_xyz=include_xyz)

        input_sys_from_mol = results_mol.get_input_system()
        input_sys_from_sys = results_sys.get_input_system()
        for input_sys in [input_sys_from_mol, input_sys_from_sys]:
            assert input_sys.has_same_atoms(system)
            assert input_sys.has_same_geometry(system)

    @pytest.mark.parametrize("include_xyz", [True, False])
    def test_get_main_molecule(self, mol, system, include_xyz):
        results_mol = self._create_results_from_molecule(mol, include_xyz=include_xyz)
        results_sys = self._create_results_from_chemical_system(system, include_xyz=include_xyz)

        main_mol_from_mol = results_mol.get_main_molecule()
        main_mol_from_sys = results_sys.get_main_molecule()

        expected_coords = mol.as_array()
        expected_lattice = np.array(mol.lattice)
        if include_xyz:
            expected_coords += self.coord_shift
            expected_lattice *= self.lattice_scaling

        for main_mol in [main_mol_from_mol, main_mol_from_sys]:
            assert np.all(main_mol.symbols == mol.symbols)
            assert np.allclose(main_mol.as_array(), expected_coords)
            assert np.allclose(main_mol.lattice, expected_lattice)
            assert main_mol.label(4) == mol.label(4)

    @pytest.mark.parametrize("include_xyz", [True, False])
    def test_get_main_system(self, mol, system, include_xyz):
        results_mol = self._create_results_from_molecule(mol, include_xyz=include_xyz)
        results_sys = self._create_results_from_chemical_system(system, include_xyz=include_xyz)

        expected_sys = system.copy()
        if include_xyz:
            expected_sys.coords += self.coord_shift
            expected_sys.lattice.vectors *= self.lattice_scaling

        main_sys_from_mol = results_mol.get_main_system()
        main_sys_from_sys = results_sys.get_main_system()
        for main_sys in [main_sys_from_mol, main_sys_from_sys]:
            assert main_sys.has_same_atoms(expected_sys)
            assert main_sys.has_same_geometry(expected_sys)

    @pytest.mark.parametrize("include_xyz", [True, False])
    def test_get_main_ase_atoms(self, mol, system, include_xyz):
        results_mol = self._create_results_from_molecule(mol, include_xyz=include_xyz)
        results_sys = self._create_results_from_chemical_system(system, include_xyz=include_xyz)

        ase_atoms_from_mol = results_mol.get_main_ase_atoms()
        ase_atoms_from_sys = results_sys.get_main_ase_atoms()

        expected_coords = mol.as_array()
        expected_lattice = np.zeros((3, 3))
        for i, vec in enumerate(mol.lattice):
            expected_lattice[i, :] = vec
        if include_xyz:
            expected_coords += self.coord_shift
            expected_lattice *= self.lattice_scaling

        for ase_atoms in [ase_atoms_from_mol, ase_atoms_from_sys]:
            assert np.all(ase_atoms.numbers == mol.numbers)
            assert np.allclose(ase_atoms.positions, expected_coords)
            assert np.allclose(ase_atoms.cell, expected_lattice)
            assert fromASE(ase_atoms).label(4) == mol.label(4)


class TestAMSWorkerResultsWithHydroxide(WorkerResultsTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "hydroxide.xyz")
        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]]
        mol.properties.charge = -1
        mol[1].properties.region = {"oxygen", "hydroxide"}
        mol[2].properties.region = {"hydrogen", "hydroxide"}
        mol[2].properties.mass = 2.0141
        mol[1].properties.adf.f = "myfrag"
        mol[2].properties.adf.f = "myfrag"
        mol.guess_bonds()
        return mol

    @pytest.fixture
    def system(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "hydroxide.xyz"))

        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]]
        mol.charge = -1
        mol.add_atoms_to_region([0, 1], "hydroxide")
        mol.add_atom_to_region(0, "oxygen")
        mol.add_atom_to_region(1, "hydrogen")
        mol.enable_atom_attributes("adf")
        mol.atoms[1].mass = 2.0141
        mol.atoms[0].adf.f = "myfrag"
        mol.atoms[1].adf.f = "myfrag"
        mol.guess_bonds()
        return mol


class TestAMSWorkerResultsWithWater(WorkerResultsTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "water.xyz")
        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0]]
        mol[1].properties.region = {"oxygen"}
        mol[2].properties.mass = 2
        mol[1].properties.forcefield.type = "o"
        mol[2].properties.forcefield.type = "h2"
        mol.guess_bonds()
        return mol

    @pytest.fixture
    def system(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "water.xyz"))

        mol.lattice = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0]]
        mol.add_atom_to_region(0, "oxygen")
        mol.enable_atom_attributes("adf")
        mol.enable_atom_attributes("forcefield")
        mol.atoms[1].mass = 2
        mol.atoms[0].forcefield.type = "o"
        mol.atoms[1].forcefield.type = "h2"
        mol.guess_bonds()
        return mol


class TestAMSWorkerResultsWithBenzene(WorkerResultsTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "benzene.xyz")
        mol.guess_bonds()
        return mol

    @pytest.fixture
    def system(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "benzene.xyz"))
        mol.guess_bonds()
        return mol


class TestAMSWorkerResultsWithChlorophyll(WorkerResultsTestBase):

    @pytest.fixture
    def mol(self, xyz_folder):
        mol = Molecule(xyz_folder / "chlorophyl1.xyz")
        return mol

    @pytest.fixture
    def system(self, xyz_folder):
        skip_if_no_scm_base()

        from scm.base import ChemicalSystem

        mol = ChemicalSystem.from_xyz(str(xyz_folder / "chlorophyl1.xyz"))
        return mol
