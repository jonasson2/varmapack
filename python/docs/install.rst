Installation
============

Varmapack is installed from PyPI. The commands below show common ways to create
or select a Python environment and install the package into it. To use
Varmapack in a Python program, begin with ``import varmapack``.

Conda environment
-----------------

.. code-block:: sh

   conda create -n varmapack python
   conda activate varmapack
   pip install varmapack

uv virtual environment
----------------------

.. code-block:: sh

   uv venv
   source .venv/bin/activate
   uv pip install varmapack

venv with pip
-------------

.. code-block:: sh

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install varmapack

Development version from GitHub
-------------------------------

To install the current development version directly from GitHub:

.. code-block:: sh

   python -m pip install "git+https://github.com/jonasson2/varmapack.git#subdirectory=python"

This builds Varmapack from source. For editable development checkouts and build
requirements, see `DEVELOPMENT.md
<https://github.com/jonasson2/varmapack/blob/main/DEVELOPMENT.md>`_.
