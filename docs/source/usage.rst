Usage
==========

Installation
---------------
Follow these instructions to clone the EUVpy github repo setup a new Python virtual environment with the required dependencies.

1. Clone the repo into the directory of your choosing

.. code-block:: console

    $ cd ~/code
    $ git clone https://github.com/DanBrandt/EUVpy.git

2. Ensure Fortran is installed and compile the codes needed to run the HEUVAC model.

.. code-block:: console

    $ sudo apt-get install gfortran
    $ git submodule update --init --recursive --remote
    $ ./install_heuvac.sh

3. Setup a Python virtual environment to run the code.

.. code-block:: console

    $ python3 -m venv ~/python_venvs/euv_venv

4. Activate virtual environment and install Python dependencies

.. code-block:: console

    $ source ~/python_venvs/euv_venv/bin/activate
    (euv_venv) $ cd EUVpy
    (euv_venv) $ pip install -r requirements.txt

5. Install EUVpy

.. code-block:: console

    (euv_venv) $ pip install EUVpy