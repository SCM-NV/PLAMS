Analysis tools: Densf, FCF, analysis
------------------------------------

Apart from main computational programs mentioned above, Amsterdam Modeling Suite offers a range of small utility tools that can be used to obtain more specific results.
These tools usually base on the prior run of one of the main programs and need the KF file produced by them as a part of the input.

From the functional point of view these tools are very similar to ADF or AMS.
Their results are stored in KF files and their input files follow the same structure of blocks, keys and values.
Because of that the same classes (|SCMJob| and |SCMResults|) are used as bases and hence preparation, running and results extraction for utility tools follow the same rules as the AMS program.

The main difference is that usually utility jobs don't need molecular coordinates as part of the input (they extract this information from previous calculation's KF file).
So no |Molecule| instance is needed and the ``molecule`` attribute of the job object is simply ignored.
Because of that :meth:`~SCMResults.get_molecule` method does not work with :class:`FCFResults`, :class:`DensfResults` etc.

Below you can find the list of dedicated job classes that are currently available.
Details about input specification for those jobs can be found in corresponding part of AMS suite documentation.


.. currentmodule:: scm.plams.interfaces.adfsuite.fcf

.. autoclass:: FCFJob(inputjob1=None, inputjob2=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type


.. currentmodule:: scm.plams.interfaces.adfsuite.densf

.. autoclass:: DensfJob(inputjob=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type


.. currentmodule:: scm.plams.interfaces.adfsuite.amsanalysis
.. autoclass:: AMSAnalysisJob(name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type

