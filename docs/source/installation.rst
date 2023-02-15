.. SHERLOCK PIPEline documentation master file, created by
   sphinx-quickstart on Thu Jul  8 08:43:51 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=============================================
Installation
=============================================

Firstly, you must know that *WATSON* is a package that lies on many different dependencies whose versions are
fixed. Hence, we strongly encourage all the users to use a *Python*
`virtual environment <https://docs.python.org/3/library/venv.html>`_ to install *WATSON* to avoid
versions collisions with your typical host installation.

Some of the *WATSON* dependencies need ``numpy`` and ``Cython``  before being installed and therefore you should install them before trying the SHERLOCK installation. Take also into account that the dependencies brought by SHERLOCK also need additional operating system tools that need to be installed separately (C++ compiler, Fortran compiler, etc). So, be aware of counting with the next OS packages installed in your environment (e.g. for Python3.8 in a Linux distribution):

.. code-block::

   build-essential
   libssl-dev
   python3.8
   python3-pip
   libbz2-dev
   libssl-dev
   libreadline-dev
   libffi-dev
   libsqlite3-dev
   tk-dev
   libpng-dev
   libfreetype6-dev
   llvm-9
   llvm-9-dev
   gfortran
   gcc
   locales
   python3-tk
   python3.8-dev

In case you are running a non-Linux distribution you will need to guess your OS packages matching the ones we mention for Linux.

Therefore, once you have got the OS tools, the *WATSON* package can be installed in several ways. The cleanest one is by running:

``python3 -m pip install -r requirements.txt``

Where the ``requirements.txt`` file should contain the same than
`requirements.txt <https://github.com/PlanetHunters/watson/blob/master/requirements.txt>`_.

An additional installation method is running the next commands:

``python3 -m pip install numpy``

``python3 -m pip install Cython``

``python3 -m pip install dearwatson``

You can also use our Docker image from `DockerHub <https://hub.docker.com/repository/docker/sherlockpipe/watson>`_
or build it from our `Dockerfile <https://github.com/PlanetHunters/watson/blob/master/docker/Dockerfile>`_.
Therefore, you can also use as a Singularity container meanwhile they support Docker.

-------------
Dependencies
-------------

All the needed dependencies should be included by your `pip` installation of *WATSON*. If you are
interested you can inspect the requirements list under
`setup.py <https://github.com/PlanetHunters/watson/blob/master/setup.py>`_.

-----------
Integration
-----------

*WATSON* integrates with several third party services. Some of them are listed below:

* TESS, Kepler and K2 databases through `Lightkurve <https://github.com/KeplerGO/lightkurve>`_, `ELEANOR <https://adina.feinste.in/eleanor/) and [LATTE](https://github.com/noraeisner/LATTE>`_.
* MAST and Vizier catalogs through `Lightkurve <https://github.com/KeplerGO/lightkurve>`_, `transitleastsquares <https://github.com/hippke/tls>`_ and `Triceratops <https://github.com/stevengiacalone/triceratops>`_

--------
Testing
--------

**WATSON** comes with a light automated tests suite that can be executed with:

``tox``

In case you want to test the entire *WATSON* functionality we encourage you to
run some (or all) the `manual examples <https://github.com/PlanetHunters/watson/tree/master/examples>`_.
If so, please read the instructions provided there to execute them.