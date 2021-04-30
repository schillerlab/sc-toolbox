.. highlight:: shell

============
Installation
============


Stable release
--------------

To install sc-toolbox, run this command in your terminal:

.. code-block:: console

    $ pip install sc-toolbox

This is the preferred method to install sc-toolbox, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for sc-toolbox can be downloaded from the `Github repo`_.
Please note that you require `poetry`_ to be installed.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/schillerlab/sc-toolbox

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/schillerlab/sc-toolbox/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ make install

Developers should also have `nox`_ and `nox-poetry`_ installed.


.. _Github repo: https://github.com/schillerlab/sc-toolbox
.. _tarball: https://github.com/schillerlab/sc-toolbox/tarball/master
.. _poetry: https://python-poetry.org/
.. _nox: https://nox.thea.codes/en/stable/
.. _nox-poetry: https://github.com/cjolowicz/nox-poetry
