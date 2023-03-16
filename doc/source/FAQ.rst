.. _FAQ:

FAQ
###


Do I always need to run PLAMS scripts through $ADFBIN/plams myscript.py?
************************************************************************

Highly advisable as this ensures you get the right /compatible python libraries loaded.
Alternatively you could import the plams packages:

::

   from scm.plams import *
   init()


Where can I find some help on setting up python workflows with PLAMS?
*********************************************************************

You can find a hands-on demonstration of how to write and execute a PLAMS workflow in our `video tutorial series <https://www.youtube.com/watch?v=1-gN6HJHseM>`_.

For more simple examples and recipes to get you started, take a look at the other parts in this PLAMS manual.

Let us know if you have issues setting up your own workflows. We would also like to hear if you would like to contribute a useful python workflow script to help out others!

Where can I find some example scripts?
**************************************

A good example of a multi-step PLAMS workflow with DFTB and BAND is given in our recent `video tutorial <https://www.youtube.com/watch?v=1-gN6HJHseM>`_.
Download: :download:`workflow.py <https://downloads.scm.com/distr/workflow.py>`

For several simple and more advanced workflows, take a look at the :ref:`Examples <examples>` section of this PLAMS manual.

Another common PLAMS application is the post-processing of calculation results. A good example of this is found in the analysis section of the `Battery Discharge tutorial <../Tutorials/MolecularDynamicsAndMonteCarlo/GCMCLiSBattery.html>`__.
Download: :download:`LiVoltageProfile.py <../../../../../userdoc/Tutorials/downloads/GCMCLiSBattery/LiVoltageProfile.py>`

