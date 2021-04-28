Usage
=====

API
---

Import the sc-tools API as follows:

.. code:: python

   import sc_toolbox.api as sct

You can then access the respective modules like:

.. code:: python

   sct.plot.cool_fancy_plot()

.. contents::
    :local:
    :backlinks: none

Plots
~~~~~

.. automodule:: sc_toolbox.api.plot
   :members:

Calculations
~~~~~~~~~~~~

.. automodule:: sc_toolbox.api.calc
   :members:

Utilities
~~~~~~~~~

.. automodule:: sc_toolbox.api.util
   :members:

CLI
---

.. click:: sc_toolbox.__main__:sc_toolbox_cli
   :prog: sc-toolbox
   :nested: full
