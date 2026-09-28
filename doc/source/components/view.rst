.. _view:

View
----

.. currentmodule:: scm.plams.tools.view

.. note::

    The |view| function is available in AMS2026+

The |view| function aims to simplify the process of visualizing molecular and periodic systems in PLAMS.
It can generate images for displaying systems in a Jupyter notebook as well as saving images to files.

Movies
~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    The |movie| and |display_movie| functions are available in AMS2027+

The |movie| function creates animations from existing image frames.
Frames can be image file paths, a directory containing image files, PIL images, matplotlib figures or matplotlib axes.
This makes it possible to combine |movie| with |view|, because |view| returns a PIL image.

For example, you can make a GIF or animated WebP from a sequence of molecules:

.. code-block:: python

    from scm.plams import movie, view

    frames = (view(mol) for mol in trajectory[::10])
    movie(frames, "trajectory.webp", fps=12)

GIF and animated WebP movies are written directly with Pillow.
Animated WebP is usually much smaller than GIF and is useful for web pages and notebooks.
MP4 output is also supported, but requires an external ``ffmpeg`` executable:

.. code-block:: python

    movie("frames/", "trajectory.mp4", fps=24)

If ``ffmpeg`` cannot be found, |movie| raises an error explaining how to provide it.
PLAMS searches for an explicitly supplied ``ffmpeg`` path, ``SCM_FFMPEG``, ``$AMSBIN/ffmpeg``, and finally ``ffmpeg`` on ``PATH``.

Use |display_movie| to display a saved movie in a Jupyter notebook:

.. code-block:: python

    from scm.plams import display_movie

    display_movie("trajectory.webp", width=600)

Image grids
~~~~~~~~~~~~~~~~~~~~~~~~~

The |plot_image_grid| function displays a dictionary of images in a matplotlib grid.
It is useful when you want to compare several images from |view| in one figure:

.. code-block:: python

    from scm.plams import plot_image_grid, view

    images = {
        "initial": view(initial_molecule),
        "final": view(final_molecule),
    }
    plot_image_grid(images, cols=2)

For a detailed worked example demonstrating the capabilities and uses of the |view| function, see the `Visualization example <../../PythonExamples/visualization/index.html>`__.
Otherwise see below for the full API specification.

API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: scm.plams.tools.view

.. autofunction:: view
.. autoclass:: ViewConfig
    :exclude-members: __eq__, __hash__, __init__, __repr__, __weakref__

.. currentmodule:: scm.plams.tools.movie

.. autofunction:: movie
.. autofunction:: display_movie

.. currentmodule:: scm.plams.tools.plot

.. autofunction:: plot_image_grid
    :noindex:
