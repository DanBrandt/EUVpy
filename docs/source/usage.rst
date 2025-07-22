Usage
==========

Installation
---------------
Follow these instructions to clone the EUVpy github repo setup a new Python virtual environment with the required dependencies.

1. Clone the repo into the directory of your choosing

.. code-block:: console

    $ cd ~/code
    $ git clone https://github.com/DanBrandt/EUVpy.git

2. Setup a Python virtual environment to run the code.

.. code-block:: console

    $ python3 -m venv ~/python_venvs/euv_venv

3. Activate virtual environment and install Python dependencies

.. code-block:: console

    $ source ~/python_venvs/euv_venv/bin/activate
    (euv_venv) $ cd EUVpy
    (euv_venv) $ pip install -r requirements.txt

4. Install EUVpy

.. code-block:: console

    (euv_venv) $ pip install EUVpy