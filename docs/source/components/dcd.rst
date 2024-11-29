DCD trajectory files
~~~~~~~~~~~~~~~~~~~~~~~~

The DCD file is the most compact trajectory file format,
which is very advantages considering the large size of MD trajectories.
It only stores the coordinates and lattice vectors,
so for visualization it needs to be combined with another file.
VMD can read a DCD file in combination with any single frame file (XYZ, PDB),
or with a PSF topology file (CHARMM).

.. autoclass :: scm.plams.trajectories.dcdfile.DCDTrajectoryFile
    :inherited-members:
    :exclude-members: _read_variable, _write_variable, _read_header, _write_header, _is_endoffile, _convert_cell, _rewind_to_first_frame, _rewind_n_frames, __iter__, __next__, __len__, __enter__, __exit__, __del__, _move_cursor_without_reading, _set_plamsmol, _set_plamsmol, _read_plamsmol, __delattr__, __dir__, __eq__, __format__, __ge__, __getattribute__, __gt__, __hash__, __init_subclass__, __le__, __lt__, __ne__, __new__, __reduce__, __reduce_ex__, __repr__, __setattr__, __sizeof__, __str__, __subclasshook__, __weakref__

