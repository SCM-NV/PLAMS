
.. _RKF_trajectory_file:

RKF trajectory files
~~~~~~~~~~~~~~~~~~~~


.. autoclass :: scm.plams.trajectories.rkffile.RKFTrajectoryFile
    :exclude-members: _read_header, _set_mddata_units, _write_header, _write_general_section, _write_molecule_section, _set_mdunits, _read_coordinates, _read_cell_data, _read_bond_data, _store_mddata_for_step, _is_endoffile, _write_history_entry, _write_bonds_in_history, _write_mdhistory_entry, _write_keydata_in_history, _get_block_info

RKF history files
+++++++++++++++++

This subsection describes the API of the |RKFHistoryFile| class, 
which can read and write the results from simulations with changing numbers of atoms.
The majority of molecular simulations explore a subspace of the canonical, micro-canonical,
or isothermal-isobaric ensembles, in which the number of atoms :math:`N` remains constant.
However, a Grand Canonical Monte Carlo simulation is one of the exceptions in which the number of atoms in the
system does change.
The |RKFTrajectoryFile| object cannot read and write the resulting simulation history,
and the derived class |RKFHistoryFile| was developed to handle these atypical trajectories.
While the methods in this class will be slower than the ones in the parent class, the API is nearly identical.
The only exception is the :meth:`write_next` method, which has an additional argument ``elements``.

.. autoclass :: scm.plams.trajectories.rkfhistoryfile.RKFHistoryFile
       :exclude-members: set_elements, _read_header, _write_header, _read_coordinates, _read_elements_for_frame, _write_version_history, _write_system_version_history_entry, _find_system_change, _check_for_chemical_system

