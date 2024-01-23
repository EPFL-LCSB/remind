A package to reconstruct metabolic interaction networks in microbial communities
using decomposed in silico Minimal exchanges

1)Use DiME (in silico Minimal Exchanges)

You should have your solver see instructions inside docker

2) latest version of pytfa to be included in your docker file
    if not change entrypoint in the docker file

3)Files are stored as .h5 binary files for this

ReMIND: Reconstruction of Microbial Interaction Networks using Decomposed in silico minimal exchanges
==========================================
This repository contains the workflow to reconstruct metabolic interaction networks in microbial communities
using decomposed in silico Minimal exchanges


More information can be found in the
The package is available for both python and matlab.
The results in the preprint are generated with the python version.

The package is developed using python 3.6 and run in Docker (20.10.6) containers.
Tested with solvers cplex (v12.8.0) and gurobi(v9.1.2) (default implemented on cplex_interface)

Recommended to be run in docker containers, with dedicated solver installed.
`Setting up the python API of cplex <https://www.ibm.com/docs/en/icos/12.8.0.0?topic=cplex-setting-up-python-api>`_  in case docker based installation is not used

Generated data used in the manuscript is available under
`data <https://github.com/EPFL-LCSB/remind/tree/master/python/remind/projects/bee_project/data>`_ subfolder
under hierarchical data format.

Requirements
------------
You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/remind.git /path/to/remind
    cd /path/to/remind
    git lfs install
    git lfs pull

Further the following pip-python packages are required (can be found in detail in requirements.txt

    - optlang
    - cobra==0.17.1
    - numpy<=1.17.4
    - pandas
    - matplotlib
    - tables
    - sklearn
    - ipython
    - jedi==0.17.2
    - tqdm
    - scipy
    - holoviews
    - matplotlib_venn
    - pytfa


Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/remind/tree/master/python/remind/docker


.. code:: bash

    cd remind/python/remind/docker
    ./build.sh
    ./run.sh

Building the docker image takes approximately 5 mins.



Setup
=====
If container-based installation is not preferred you can also install this module from source using ``pip``:
*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    git clone https://github.com/EPFL-LCSB/remind.git /path/to/remind/python
    pip3 install -e /path/to/remind

The installation process should not exceed a minute if the requirements are installed. If they are not, it might take longer as the installer installs them first.

