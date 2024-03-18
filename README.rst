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

The package is developed using python 3.6 and run in Docker (20.10.21) containers.
Tested with solvers cplex (v12.8.0) and gurobi(v9.1.2)

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

Quick start
===========
As a demo following examples can be run for the 2-member honeybee gut community after building the environment or insid ethe docker
This tutorial aims to show the step to use our framework. Can be adapted to any community with the desired extracellular environment.

.. code:: bash
        cd /
        cd remind/projects/tutorial/

First get the DiMEs for both members by running the following bash script. Number of alternatives are limited to 10
for tutorial purposes can be changed inside the `get_dimes_tutorial.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/get_dimes_tutorial.py>`_ script.
script get_dimes_tutorial.py by modifying the max_alternative.
    .. code:: bash
    ./bash_tutorial_dimes.sh

After generating the DiMEs merge the DiMEs and build the community model with the following script inside Ipython.
.. code-block:: python
    run build_community_model_from_dimes_tutorial.py

The next step is to use the built community model and reconstruct the interaction networks with a user defined objective function
via the ILP formulation. For this you can refer to the `run_ilp_tutorial_community_model.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/run_ilp_tutorial_community_model.py>`_ script.
for various objective functions. To run for the indicated objective functions run the following bash script.
    .. code:: bash
    ./bash_tutorial_ilp.sh

After running the ILP for various objective functions you can analyse the data inside Ipython:
.. code-block:: python
    run analysis_ilp_solutions_tutorial.py
    #check the alternative cooperation patterns
     print(frame_int_coop.pos_int)

