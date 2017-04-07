You must be in the right python enviornment (2.7, certain packages
installed, clear of mind, pure of soul, etc.)

To accomplish this, open your anaconda terminal (download
[anaconda](https://www.continuum.io/downloads) or
[miniconda](https://conda.io/miniconda.html)), navigate to this
folder, and type the following command:

```
conda env create python=2.7
```

This will create a conda environment, downloading things as neccesary,
which is suitable to run CATS.

It might take a few minutes.  Go have a coffee.

When it is done, you must activate the new python envionrment each
time, before running cats.  Basically, it tells the anaconda terminal
which version of python you want to run.

To do this, type the following command:

```
activate cats_python
```

Then, you should be all setup to go!  Cats is run using the command:

```
python cats_v1.py
```

Happy tracking of combinatorial states :)

 - Dan
