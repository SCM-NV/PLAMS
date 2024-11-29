from scm.plams.recipes.adfcosmorsconformers import ADFCOSMORSConfJob, ADFCOSMORSConfFilter
from scm.plams import from_smiles, Settings
from scm.conformers import ConformersJob

mol = from_smiles("CC(=O)O")

conf_sett = Settings()
conf_sett.input.AMS.Generator.RDKit
conf_sett.input.AMS.Generator.RDKit.InitialNConformers = 50
conf_job = ConformersJob(name="conformers_uff", molecule=mol, settings=conf_sett)

dftb_sett = Settings()
dftb_sett.input.AMS.Task = "Optimize"
dftb_sett.input.DFTB

# ADFCOSMORSConfFilter(max number of conformers, max energy range)
fil1 = ADFCOSMORSConfFilter(20, 22)  # applied to UFF
fil2 = ADFCOSMORSConfFilter(10, 12)  # applied to DFTB
fil3 = ADFCOSMORSConfFilter(5, 7)  # applied to ADF gas phase

a = ADFCOSMORSConfJob(
    mol,
    conf_gen=conf_job,
    first_filter=fil1,
    additional=[(dftb_sett, fil2)],
    final_filter=fil3,
    coskf_name="acetic_acid",
    coskf_dir="test_coskfs",
)
a.run()
