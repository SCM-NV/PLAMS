
.. _trajectory_class:

Trajectory class
~~~~~~~~~~~~~~~~~~~~

This subsection describes the API of the |Trajectory| class.
While the |RKFTrajectoryFile| clearly represents a file, the |Trajectory| is set up to represent the trajectory itself.
It behave like a list of |Molecule| objects,
while keeping the memory requirements to a minimum.
A |Trajectory| object can be associated with a single RKF file, or with multiple RKF files.
In the latter case, the frames will be concatenated.

.. autoclass :: scm.plams.trajectories.trajectory.Trajectory
    :exclude-members: _get_filenum_and_stepnum, _get_index, __init__, __len__, __iter__, __getitem__, __weakref__


