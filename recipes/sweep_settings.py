from pathlib import Path
from typing import List
from scm.plams import Settings


def make_sweep(s_default: Settings, parameter_interface: Settings):
    """make a sweep over parameter. Note each settings key is independent to each other:

    parameter_interface = Settings()
    parameter_interface.input.seed = [1, 2]
    parameter_interface.input.temperature = [300]
    output = make_sweep(s_default = Settings(), parameter_interface=parameter_interface)

    # it results in:
    s = Settings()
    s.input.seed = 1
    print(output[0] == s)
    s = Settings()
    s.input.seed = 2
    print(output[1] == s)
    s = Settings()
    s.input.temperature = 300
    print(output[2] == s)

    :param s_default: default setting on which make the sweep
    :type s_default: Settings
    :param parameter_interface: the settings with each value should be a list
    :type parameter_interface: Settings
    :raises TypeError: if the settings values are not list a TypeError is raised
    :return: list of settings
    :rtype: List[Settings]
    """
    settings_collection = []
    parameter_sweep = parameter_interface.flatten(flatten_list=False).as_dict()
    for path_key, list_values in parameter_sweep.items():
        if not isinstance(list_values, list):
            raise TypeError(f"For {path_key=} the {type(list_values)=} but must be List")
        for values_i in list_values:
            sett = Settings()
            sett.set_nested(path_key, values_i)
            sett += s_default.copy()
            settings_collection.append(sett)
    return settings_collection


def grid_search(s_defaults: List[Settings], parameter_interface: Settings):
    """
    prepares a grid search among some parameter list given by the parameter_interface
    s_defaults = [Settings()]
    parameter_interface = Settings()
    parameter_interface.input.seed = [1, 2]
    parameter_interface.input.temperature = [300, 700]

    it will returns 4 settings:
        s_defaults.input.seed s_defaults.input.temperature
        1                       300
        2                       300
        1                       700
        2                       700

    """
    for path_key, list_values in parameter_interface.flatten(flatten_list=False).as_dict().items():
        param_interface_single_path = Settings()
        param_interface_single_path.set_nested(path_key, list_values)
        new_defaults = []
        for s_default_i in s_defaults:
            s_defaults_out = make_sweep(s_default=s_default_i.copy(), parameter_interface=param_interface_single_path)
            new_defaults.extend(s_defaults_out)
        s_defaults = new_defaults
    return s_defaults