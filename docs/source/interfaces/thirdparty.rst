Other programs
-------------------------

In this chapter we present a list of PLAMS interfaces to other programs and packages.
A majority of what follows below are interfaces to so called *external binaries* -- computational chemistry tools that come in a form of executable programs that read an input file, perform calculations and produce various output files.
From PLAMS perspective each such interface is just a subclass of |SingleJob| that implements a way of producing a working runscript (|get_runscript|) and a valid input file (|get_input|) for a particular binary based on the contents of job settings.
Usually such a "specialized" |SingleJob| subclass comes together with a corresponding specialized |Results| subclass providing methods for accessing the data produced by the binary.
Some of these methods are just simple convenience shortcuts (like ``get_energy()`` or ``get_main_molecule()``), others provide access to files in whatever formats a particular binary produces (XML, binary files).

Interfaces described below are divided into interfaces to programs and tools that are included in Amsterdam Modeling Suite and interfaces to third party computational chemistry packages (usually contributed by other PLAMS users).
The last chapter presents a bit different kind of interfaces, so called *molecule interfaces*.
They offer a way of using PLAMS |Molecule| class with other libraries capable of manipulating molecular coordinates.

.. toctree::

    cp2k
    crystal
    dftbplus
    dirac
    orca
    raspa
    vasp
