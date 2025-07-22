Documentation
==============

The EUVpy documentation is built using `Sphinx <http://sphinx-doc.org/>`_. If you edit any of the .rst files in the `docs/source` directory, you will need to rebuild the documentation.
To do this you will need to install Sphinx and the required themes in your Python virtual environment.

.. code-block:: console

    $ source venv/bin/activate
    $ pip install sphinx

Now you can rebuild the documentation by changing into the EUVpy `docs` directory and running:

.. code-block:: console

    $ make clean html
    $ make html

Now if you open the docs in your browser you should see your changes.  One way to open the docs is to navigate to the `docs/build/html` directory and open the `index.html` file in your browser.