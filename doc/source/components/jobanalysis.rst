.. _job_analysis:

Job Analysis
-------------

.. currentmodule:: scm.plams.tools.job_analysis

.. note::

    The |JobAnalysis| class is available in AMS2025+

The |JobAnalysis| class is a tool that aims to simplify the process of analyzing the status, inputs and outputs of multiple jobs.
It helps to create analysis tables, with jobs as rows and "analysis fields" as columns.
An analysis field is simply a definition for how to extract a value from each job.
These tables can then be easily visualized in a Jupyter notebook or script, or exported for use with other analysis packages like `pandas <https://pandas.pydata.org/docs/index.html>`_.

For a worked example demonstrating the capabilities and uses of the |JobAnalysis| class, see :ref:`JobAnalysisExample`.

Adding Jobs
~~~~~~~~~~~

Jobs can be added to a |JobAnalysis| in two ways, either passed directly on initialization or using :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_job`,
or alternatively loaded from paths using :meth:`~scm.plams.tools.job_analysis.JobAnalysis.load_job`.

When loading jobs from a path, an attempt is first made to load from a ``.dill`` file.
If this fails, the tool will attempt to use a series of loaders to create the job.
By default this will be :meth:`~scm.plams.interfaces.adfsuite.ams.AMSJob.load_external`, but a custom set of loaders can be provided for other job types.

For example, the following snippet will add two jobs to the |JobAnalysis| and then load a further job using the ``ParAMSJob.load_external`` method:

.. code-block:: python

    ja = (JobAnalysis(jobs=[job1, job2])
          .load_job("path/to/job3", loaders=[ParAMSJob.load_external]))

This also illustrates one of the design features of |JobAnalysis|, each operation returns an updated copy of the |JobAnalysis|, allowing the use of fluent syntax to chain methods.

Adding Analysis Fields
~~~~~~~~~~~~~~~~~~~~~~

Analysis fields can be added to a |JobAnalysis| in a number of ways.
Firstly, there are a small number of predefined "standard" analysis fields which are common across jobs e.g. ``Name``, ``Path`` etc.
These can be added to the analysis via the methods, such as :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_standard_field` / :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_standard_fields`.
To add these standard fields, the key(s) of the relevant fields must be supplied. A full list of these is available in the tooltip.

Secondly, custom analysis fields can be added using the methods :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_field` or :meth:`~scm.plams.tools.job_analysis.JobAnalysis.set_field`.
When doing this, a unique identifier for the field (the key) must be provided, along with a function defining how to extract a value for the field from a job.
In addition, optional arguments can be provided to set the displayed name for the field, and its formatting when returned in a table.

Finally, additional methods are present to facilitate adding fields for values from the job ``settings``.
These are :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_settings_field` and :meth:`~scm.plams.tools.job_analysis.JobAnalysis.add_settings_input_fields`.
The former is used to create an analysis field given a specific nested settings key, the latter automatically creates fields for all settings keys under the ``settings.input``, which includes AMS input settings.
The key for these settings fields will be the concatenated settings key in Pascal case e.g. ``("input", "ams", "task")`` will have field key ``InputAmsTask``.

For example, the following snippet will add the ``Formula`` field to the |JobAnalysis|, followed by a custom field displaying the number of atoms, and finally the settings input fields:

.. code-block:: python

    ja = (ja
     .add_standard_field("Formula")
     .add_field("NAtoms", lambda j: len(j.molecule))
     .add_settings_input_fields())


Modifying Analysis
~~~~~~~~~~~~~~~~~~

|JobAnalysis| supports some modification and data manipulation methods.
These allow jobs and results to be interrogated and more easily visualized.

Fields can be filtered using :meth:`~scm.plams.tools.job_analysis.JobAnalysis.filter_fields`.
This accepts a predicate as an argument, which accepts all the values for a field.
Fields are only retained for which the predicate evaluates to ``True``.
There are some pre-configured predicates, such as :meth:`~scm.plams.tools.job_analysis.JobAnalysis.remove_empty_fields` and :meth:`~scm.plams.tools.job_analysis.JobAnalysis.remove_uniform_fields`.
These can be useful for removing noise from the analysis table, for fields which only have empty/the same values.

Jobs can also be filtered using :meth:`~scm.plams.tools.job_analysis.JobAnalysis.filter_jobs`.
This accepts a predicate as an argument, which accepts all the field keys and values for a job.
Again, jobs are only retained for which the predicate evaluates to ``True``.
This can be useful for removing jobs which are not of interest, for example those with a given status.

For example, the following snippet will retain only jobs which did not succeed:

.. code-block:: python

    ja = (ja
     .filter_field(lambda vals: all([v is not None for v in vals]))
     .filter_jobs(lambda data: not data["OK"]))


Fields can also be sorted, renamed and formatted to aid with presentation.
These changes will only be reflected in the generated tables, not when calling :meth:`~scm.plams.tools.job_analysis.JobAnalysis.get_analysis`.

For example:

.. code-block:: python

    ja = (ja
     .rename_field("InputAmsTask", "Task")
     .format_field("Energy", ".4f")
     .sort_fields(["Name", "Formula", "Energy"])
     .sort_jobs(field_keys=["Energy"]))


Extracting Analysis
~~~~~~~~~~~~~~~~~~~

There are various ways to extract data from |JobAnalysis|.
The simplest way to get a visual representation is to call :meth:`~scm.plams.tools.job_analysis.JobAnalysis.get_table` (or :meth:`~scm.plams.tools.job_analysis.JobAnalysis.display_table` if running in a notebook).
This generates a table in either ``markdown``, ``html`` or ``rst`` formats.

Alternatively, data can be retrieved in code by calling :meth:`~scm.plams.tools.job_analysis.JobAnalysis.get_analysis`.
This gets a pure Python dictionary of the data, where the analysis field keys are the dictionary keys.
The dictionary values are lists of values, with the elements as the field value for each job.

For more complex or involved analysis, the best approach is to export the data to a pandas data frame.
Pandas is a fast and powerful data analysis tool, which can perform complex manipulations.
It can be installed via ``amspackages``.

As a final option, data can also be saved to a CSV file using :meth:`~scm.plams.tools.job_analysis.JobAnalysis.to_csv_file`.


API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: JobAnalysis
    :exclude-members: __weakref__, _Field, _get_field_analysis, _is_empty_value, _get_job_settings, _add_job, _load_job, _add_field, _set_field, _remove_field, _add_standard_field, _read_rkf, _repr_html_, __dir__, _
