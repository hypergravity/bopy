.. bopy_doc documentation master file, created by
   sphinx-quickstart on Sat Jan  2 20:46:37 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to bopy_doc's documentation!
====================================

Contents:

.. toctree::
   :maxdepth: 2

   intro
   tutorial
   TODOs


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. code-block:: python
   :emphasize-lines: 3,5

   def some_function():
       interesting = False
       print 'This line is highlighted.'
       print 'This one is not...'
       print '...but this one is.'