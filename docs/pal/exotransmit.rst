***************************************
Exo-Transmit (`ExoCTK.pal.exotransmit`)
***************************************

ExoCTK provides python wrappers and useful helper functions for
Exo-Transmit, a software package to calculate exoplanet transmission
spectra for planets of varied composition.

From the command line
---------------------

Exo-Transmit requires a working directory with opacity, equation of
state and T-P profile tables to run.  This directory can be downloaded from the web from
the command line with::

    collect_exotransmit_data

This will download opacity, equation of state, and T-P profile tables into ``Opac/``, ``EOS/``, and ``T_P/``
directories and creates a ``Spectra/`` directory for output.  It also downloads the ``userInput.in``,
``selectChem.in``, and ``otherInput.in`` files which set all of the necessary parameters for running the code.

In general, the user should only need to edit parameters in the ``userInput.in`` file.  When first downloaded some
default parameters are set.  With these files downloaded one simply needs to run::

    exotransmit

to create a transmission spectrum, which will be saved in ``Spectra/default.txt``.

From Python
-----------

Alternatively, one can run the code and create the necessary input files from within Python.
This can be useful for batch creation of transmission spectra.
With the directories described above already downloaded one can use::

    from ExoCTK.pal import exotranmsit
    exotransmit.exotransmit(output_file='from_python.txt')

This will overwrite the current ``userInput.in`` file changing the line specifying the output file and run Exo-Transmit.
The resulting spectrum will be saved in ``Spectra/from_python.txt``.
For a full list of the configurable parameters see `~ExoCTK.pal.exotransmit.create_user_input`.

.. automodapi:: ExoCTK.pal.exotransmit
    :skip: urljoin
