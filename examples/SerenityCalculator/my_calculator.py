from scm.plams.interfaces.thirdparty.serenity import SerenitySettings
from scm.plams.interfaces.thirdparty.serenity_calculator import SerenityCalculator


def get_calculator() -> SerenityCalculator:
    s = SerenitySettings()
    s.input.system.A.charge = "0"
    s.input.system.A.spin = "0"
    s.input.system.A.method = "HF"
    s.input.system.A.basis.label = "6-31GS"
    s.input.task.SCF.act = "A"

    calc = SerenityCalculator(settings=s)
    return calc
