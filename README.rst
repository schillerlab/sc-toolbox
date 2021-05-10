.. image:: https://user-images.githubusercontent.com/21954664/116578141-65a85180-a911-11eb-9c33-925a2ec600c6.png
    :target: https://github.com/schillerlab/sc-toolbox
    :alt: sc-toolbox logo
    :align: center
    :width: 200px


sc-toolbox
==========

|PyPI| |Downloads| |Python Version| |License| |Read the Docs| |Build| |Tests| |Codecov| |pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/sc-toolbox.svg
   :target: https://pypi.org/project/sc-toolbox/
   :alt: PyPI
.. |Downloads| image:: https://pepy.tech/badge/sc-toolbox
    :target: https://pepy.tech/badge/sc-toolbox
    :alt: PyPI Downloads
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/sc-toolbox
   :target: https://pypi.org/project/sc-toolbox
   :alt: Python Version
.. |License| image:: https://img.shields.io/github/license/schillerlab/sc-toolbox
   :target: https://opensource.org/licenses/MIT
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/sc-toolbox/latest.svg?label=Read%20the%20Docs
   :target: https://sc-toolbox.readthedocs.io/
   :alt: Read the documentation at https://sc-toolbox.readthedocs.io/
.. |Build| image:: https://github.com/schillerlab/sc-toolbox/workflows/Build%20sc-toolbox%20Package/badge.svg
   :target: https://github.com/schillerlab/sc-toolbox/actions?workflow=Package
   :alt: Build Package Status
.. |Tests| image:: https://github.com/schillerlab/sc-toolbox/workflows/Run%20sc-toolbox%20Tests/badge.svg
   :target: https://github.com/schillerlab/sc-toolbox/actions?workflow=Tests
   :alt: Run Tests Status
.. |Codecov| image:: https://codecov.io/gh/schillerlab/sc-toolbox/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/schillerlab/sc-toolbox
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

.. warning::
    This package is still under heavy development. It is primarily designed for in-house analyses at the `Theislab <https://github.com/theislab>`_
    and `Schillerlab <https://github.com/schillerlab>`_. Don't yet expect it to be well tested or documented.
    However, contributions are highly welcome and we will provide guidance if required.


Features
--------

* Create minimal best-practice containers for single-cell data analysis with Scanpy
* API for advanced Scanpy based plots and analyses

.. figure:: https://user-images.githubusercontent.com/21954664/116225631-5fb84200-a752-11eb-9489-16571428918f.png
   :alt: Preview plot

.. figure:: https://user-images.githubusercontent.com/21954664/116225765-824a5b00-a752-11eb-8cbf-c14ebd9ac030.png
   :alt: Preview plot 2

.. figure:: https://user-images.githubusercontent.com/21954664/116226005-c5a4c980-a752-11eb-9846-8dc72315d373.png
   :alt: Preview plot 3

Installation
------------

You can install *sc-toolbox* via pip_ from PyPI_:

.. code:: console

   $ pip install sc-toolbox

Usage
-----

.. code:: python

   import sc_toolbox.api as sct

Please see the `Usage documentation <Usage_>`_.

Credits
-------

This package was created with cookietemple_ using cookiecutter_ based on Hypermodern_Python_Cookiecutter_.
Most scripts were developed by `Meshal Ansari <https://github.com/mesh09/>`_ and the package was designed by `Lukas Heumos <https://github.com/zethson>`_.

.. _cookietemple: https://cookietemple.com
.. _cookiecutter: https://github.com/audreyr/cookiecutter
.. _MIT: http://opensource.org/licenses/MIT
.. _PyPI: https://pypi.org/
.. _Hypermodern_Python_Cookiecutter: https://github.com/cjolowicz/cookiecutter-hypermodern-python
.. _pip: https://pip.pypa.io/
.. _Usage: https://sc-toolbox.readthedocs.io/en/latest/usage.html
.. _API: https://sc-toolbox.readthedocs.io/en/latest/api.html
