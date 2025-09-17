.. _public-functions:

Public functions
-------------------------

.. currentmodule:: scm.plams.core.functions

This chapter gathers information about public functions that can be used in PLAMS scripts.

.. autofunction:: init
.. autofunction:: finish
.. autofunction:: load
.. autofunction:: load_all
.. autofunction:: read_molecules
.. autofunction:: read_all_molecules_in_xyz_file

.. _logging:

Logging
~~~~~~~~~~~~~~~~~~~~~~~~~

PLAMS features a simple logging mechanism.
All important actions happening in functions and methods register their activity using log messages.
These massages can be printed to the standard output and/or saved to the logfile located in the main working folder.

Every log message has its "verbosity" defined as an integer number: the higher the number, the more detailed and descriptive the message is.
In other words, it is a measure of importance of the message.
Important events (like "job started", "job finished", "something went wrong") should have low verbosity, whereas less crucial ones (for example "pickling of job X successful") -- higher verbosity.
The purpose of that is to allow the user to choose how verbose the whole logfile is.
Each log output (either file or stdout) has an integer number associated with it, defining which messages are printed to this channel (for example, if this number is 3, all messages with verbosity 3 or less are printed).
That way using a smaller number results in the log being short and containing only the most relevant information while larger numbers produce longer and more detailed log messages.

The behavior of the logging mechanism is adjusted by ``config.log`` settings branch with the following keys:

*   ``file`` (integer) -- verbosity of logs printed to the ``logfile`` in the main working folder.
*   ``stdout`` (integer) -- verbosity of logs printed to the standard output.
*   ``time`` (boolean) -- print time of each log event.
*   ``date`` (boolean) -- print date of each log event.

Log messages used within the PLAMS code use four different levels of verbosity:

*   **1**: important
*   **3**: normal
*   **5**: verbose
*   **7**: debug

Even levels are left empty for the user.
For example, if you find level 5 too verbose and still want to be able to switch on and off log messages of your own code, you can log them with verbosity 4.

.. note::

    Your own code can (and should) contain some |log| calls.
    They are very important for debugging purposes.

.. autofunction:: log



.. _binding-decorators:

Binding decorators
~~~~~~~~~~~~~~~~~~

Sometimes one wants to expand the functionality of a class by adding a new method or modifying an existing one.
It can be done in a few different ways:

*   One can go directly to the source code defining the class and modify it there before running a script.
    Such a change is global -- it affects all the future scripts, so in most cases it is not a good thing (for defining |prerun|, for example).
*   Creating a subclass with new or modified method definitions is usually the best solution.
    It can be done directly in your script before the work is done. The newly defined class can be then used instead of the old one.
    However, this solution fails in some rare cases when a method needs to differ for different instances or when it needs to be changed during the runtime of the script.
*   PLAMS binding decorator |add_to_instance| can be used.

Binding decorators associate methods to existing class instances without having to define a subclass.
Such changes are visible only inside the script in which they are used.

The usage of the PLAMS |add_to_instance| binding decorator is straightforward.
You simply define a regular function somewhere inside your script and decorate it with |add_to_instance|.
The function needs to have a valid method syntax, so it should have ``self`` as the first argument and use it to reference the class instance.

.. autofunction:: add_to_instance

.. technical::

    The above decorator is in fact a decorator factory that, given an object (class or instance), produces a decorator that binds the function as a method of that object.
    The decorator is for adding instance methods only, it cannot be used for static or class methods.
