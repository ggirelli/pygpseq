pyGPSeq v2.0.1
=======================

A Python package that provides tools to analyze images of GPSeq samples.

Sample scripts are available showing the code for single/batch runs.

Read the [documentation](https://github.com/ggirelli/gpseq-img-py/wiki) for more details.

Installation
-------------

To **install**, run the following:

```
git clone http://github.com/ggirelli/gpseq-img-py
cd gpseq-img-py
sudo -H pip3 install -e .
```

To **uninstall** run the following from within the repository folder:

```
sudo -H python3 setup.py develop --uninstall
for f in $(ls bin); do whereis $f | cut -f 2 -d " " | xargs sudo -H rm; done
```

To **update**, first uninstall, and then run the following from within the repository folder.

```
git pull
sudo -H pip3 install -e .
```

License
---

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```