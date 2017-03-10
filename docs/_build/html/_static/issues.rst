Known issues
============

1. Font-cache is built at every run.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If everytime that pyGPSeq is imported, or everytime you run ``pygpseq.Main.run()`` you get the following (or similar) warning:

``/usr/local/lib/python2.7/dist-packages/matplotlib/font_manager.py:273:``
``UserWarning: Matplotlib is building the font cache using fc-list. This may take a moment.``

Then, try removing the ``~/.cache/matplotlib/`` folder. This should fix the issue.
