ReMIND: Reconstruction of Microbial Interaction Networks using Decomposed in silico minimal exchanges
==========================================
This repository contains the workflow to reconstruct metabolic interaction networks in microbial communities
using decomposed in silico Minimal exchanges


More information can be found in the
The package is available for both python and matlab.
The results in the preprint are generated with the python version. The details on the matlab version and setup can be found in
`this file <https://github.com/EPFL-LCSB/remind/blob/master/matlab/remind_matlab.txt>`_
The package is developed using python 3.6 and run in Docker (20.10.21) containers.
Tested with solvers cplex (v12.8.0) and gurobi(v9.1.2) and tested on both Linux and macOS operating systems.

Recommended to be run in docker containers, with dedicated solver installed.
`Installing cplex <https://www.ibm.com/support/pages/cplex-optimization-studio-v128>`_  and
`Setting up the python API of cplex <https://www.ibm.com/docs/en/icos/22.1.0?topic=cplex-setting-up-python-api>`_  in case docker based installation is not used

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

    - optlang==1.6.0
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
up. Please clone the `pytfa repository <https://github.com/EPFL-LCSB/pytfa>`_  under the docker folder, or change the entrypoint and pip install it.

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
    pip install -e /path/to/remind

The installation process should not exceed a minute if the requirements are installed. If they are not, it might take longer as the installer installs them first.
Please make sure you are installing it on your virtual environment activated with python 3.6. In case you run into an issue with swiglpk package (which is required by optlang) possible workaround is
to install it via conda-forge
.. code:: bash
    conda install conda-forge::swiglpk
Quick start - Demo
===========
As a demo following examples can be run for the 2-member honeybee gut community after building the environment or inside the docker
This tutorial aims to show the step to use our framework. Can be adapted to any community with the desired extracellular environment.
As mentioned above in these scripts most data is saved as .h5 binary files. For this you will need hdf files downloaded if you are
running inside docker. You can find the instructions in `instructions_hdf5.txt <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/docker/HDF5Files/instructions_hdf5.txt>`_ file.
Or change the storing in the scripts `get_dimes_tutorial.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/get_dimes_tutorial.py>`_
`build_community_model_from_dimes_tutorial.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/build_community_model_from_dimes_tutorial.py>`_
and
`run_ilp_tutorial_community_model.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/run_ilp_tutorial_community_model.py>`_
from "to_hdf" to another format  (e.g. csv) "to_csv".


.. code:: bash


    cd /
    cd remind/projects/tutorial/


First get the DiMEs for both members by running the following bash script. Number of alternatives are limited to 10
for tutorial purposes can be changed inside the `get_dimes_tutorial.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/get_dimes_tutorial.py>`_ script.
script get_dimes_tutorial.py by modifying the max_alternative. The time to generate alternatives for each yield regime should take around 1-2 minutes for 10 alternatives. In total running the bash file for 10 biomass yield regimes for two models should take around 20 minutes. You can comment in the bash file to generate only for certain biomass yields.


.. code:: bash


    ./bash_tutorial_dimes.sh

After generating the DiMEs merge the DiMEs and build the community model and save it with the following script inside Ipython.
Generating the community model should take around 2-5 minutes.


.. code-block:: python

    ipython
    run build_community_model_from_dimes_tutorial.py

The next step is to use the built community model and reconstruct the interaction networks with a user defined objective function
via the ILP formulation. For this you can refer to the `run_ilp_tutorial_community_model.py <https://github.com/EPFL-LCSB/remind/blob/master/python/remind/projects/tutorial/run_ilp_tutorial_community_model.py>`_ script.
for various objective functions. To run for the indicated objective functions run the following bash script, which runs for three different objective functions, and generates alternative patters based on different integer cuts, which is commented in the bash file.
For each run the simulation takes around 2-10 minutes depending on the objective funciton (e.g. depending on the number of alternatives to exhaust)


.. code:: bash


    ./bash_tutorial_ilp.sh

After running the ILP for various objective functions you can analyse the data inside Ipython, this step should take less than a minute.


.. code-block:: python

    ipython
    run analysis_ilp_solutions_tutorial.py
    #check the alternative maximal cooperation patterns
    print(frame_int_coop.pos_int)
    #check the alternative minimal competition patterns
    print(frame_int_comp.neg_int)
    #check the alternative minimal environments
    print(frame_int_abiot.abiotic)



Reproduction of Results
=======
To then generate the figures in the manuscript you can check the scripts inside the `figures <https://github.com/EPFL-LCSB/remind/tree/master/python/remind/projects/bee_project/figures>`_ subfolder.


License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/remind/blob/master/LICENSE.txt>`_ file for more details.

