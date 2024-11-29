XYZ trajectory files
~~~~~~~~~~~~~~~~~~~~~~~~


.. autoclass :: scm.plams.trajectories.xyzfile.XYZTrajectoryFile
    :inherited-members:
    :exclude-members: _read_header, _write_header, _read_coordinates, _is_endoffile, _write_moldata, _rewind_to_first_frame, _rewind_n_frames, __iter__, __next__, __len__, __enter__, __exit__, __del__, _move_cursor_without_reading, _set_plamsmol, _set_plamsmol, _read_plamsmol, __delattr__, __dir__, __eq__, __format__, __ge__, __getattribute__, __gt__, __hash__, __init_subclass__, __le__, __lt__, __ne__, __new__, __reduce__, __reduce_ex__, __repr__, __setattr__, __sizeof__, __str__, __subclasshook__, __weakref__, _convert_cell

XYZ history files
+++++++++++++++++

This subsection describes the API of the |XYZHistoryFile| class, 
which can read and write the results from simulations with changing numbers of atoms.
The majority of molecular simulations explore a subspace of the canonical, micro-canonical,
or isothermal-isobaric ensembles, in which the number of atoms :math:`N` remains constant.
However, a Grand Canonical Monte Carlo simulation is one of the exceptions in which the number of atoms in the
system does change.
The |XYZTrajectoryFile| object cannot read and write the resulting simulation history,
and the derived class |XYZHistoryFile| was developed to handle these atypical trajectories.
While the methods in this class will be slower than the ones in the parent class, the API is nearly identical.
The only exception is the :meth:`write_next` method, which has an additional argument ``elements``.

.. autoclass :: scm.plams.trajectories.xyzhistoryfile.XYZHistoryFile
       :exclude-members: _is_endoffile, _read_coordinates 

