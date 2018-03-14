pyGPSeq v2.0.0
=======================

A Python package that provides tools to analyze images of GPSeq samples.

Sample scripts are available showing the code for single/batch runs.

Read the [documentation](https://github.com/ggirelli/gpseq-img-py/wiki) for more details.

Installation
-------------

```
git clone http://github.com/ggirelli/gpseq-seq-py
cd gpseq-seq-py
sudo -H pip3 install -e .
```

To update, run the following from within the repository folder.

```
sudo -H python3 setup.py develop --uninstall
for f in $(ls bin); do whereis $f | cut -f 2 -d " " | xargs sudo -H rm; done
git pull
sudo -H pip3 install -e .
```

To uninstall run the following from within the repository folder:

```
sudo -H python3 setup.py develop --uninstall
for f in $(ls bin); do whereis $f | cut -f 2 -d " " | xargs sudo -H rm; done
```

License
---

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```