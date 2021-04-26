sc-toolbox
===========================

|PyPI| |Python Version| |License| |Read the Docs| |Build| |Tests| |Codecov| |pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/sc-toolbox.svg
   :target: https://pypi.org/project/sc-toolbox/
   :alt: PyPI
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
    This package is still under heavy development. It is primarily designed for in-house analyes at the `Theislab <https://github.com/theislab>`_
    and `Schillerlab <https://github.com/schillerlab>`_. Don't yet expect it to be well tested or documented.
    However, contributions are highly welcome and we will provide guidance if required.


Features
--------

* Create minimal best-practice containers for single-cell data analysis with Scanpy
* API for advanced Scanpy based plots and analyses

![image](https://user-images.githubusercontent.com/21954664/116152289-01ee1080-a6e6-11eb-9aeb-80f6c7c98e07.png)


Installation
------------

You can install *sc-toolbox* via pip_ from PyPI_:

.. code:: console

   $ pip install sc-toolbox


Usage
-----

Please see the `Command-line Reference <Usage_>`_ for details.


Contributing
------------

Contributions are very welcome. To learn more, see the `Contributor Guide`_.


Credits
-------

This package was created with cookietemple_ using Cookiecutter_ based on Hypermodern_Python_Cookiecutter_.

.. _cookietemple: https://cookietemple.com
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _MIT: http://opensource.org/licenses/MIT
.. _PyPI: https://pypi.org/
.. _Hypermodern_Python_Cookiecutter: https://github.com/cjolowicz/cookiecutter-hypermodern-python
.. _pip: https://pip.pypa.io/
.. _Contributor Guide: CONTRIBUTING.rst
.. _Usage: https://sc-toolbox.readthedocs.io/en/latest/usage.html
