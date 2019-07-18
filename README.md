Overview of SAEGUS
===================================

SAEGUS (Simulation and Evolution of Genomes Under Selection) is a Python package which is built around the forward-time population genetics simulator [simuPOP](http://simupop.sourceforge.net/Main/HomePage). SAEGUS adds functions for modeling quantitative traits, customized breeding schemes and handling data. SAEGUS only supports Python 3. This readme will show you how to install SAEGUS. A [user guide](https://saegus.readthedocs.io/en/latest/index.html) is provided separately.

Credits
===================================

SAEGUS was developed and coded by John J. Dougherty III under the guidance of Dr. Randall J. Wisser. This work was supported by the Agriculture and Food Research Initiative Grant No. 2011-67003-30342 from the USDA-NIFA; the National Science Foundation Grant No. 1144726, Systems Biology of Cells in Engineered Environments IGERT program; and University of Delaware's Center for Bioinformatics and Computational Biology graduate program.

Installing the Anaconda Distribution
================================================

[Anaconda](https://www.anaconda.com/distribution/) <!-- https://conda.io/docs/user-guide/install/index.html) --> is a package manager and a Python distribution which ships with popular Python packages used in scientific computing. The Anaconda distribution must be installed in order to proceed with this readme.

Installing simuPOP and h5py
====================================================

[simuPOP](http://simupop.sourceforge.net/) can be installed using the conda package manager. Once you have installed Anaconda, open a terminal and run:

```bash
$ conda install -c conda-forge simuPOP
```

Having toruble? See source documentation on installation: [simuPOP Installation](https://github.com/BoPeng/simuPOP)<!-- (https://anaconda.org/bpeng/simupop) -->

[h5py](https://www.h5py.org/) is a Python package for hdf5 files. I am unsure if it ships by default with the Anaconda distribution. To be sure, you can run the install command regardless.

```bash
$ conda install -c anaconda h5py
```
Having toruble? See source documentation on installation: [h5py Installation](http://docs.h5py.org/en/stable/build.html)


Setting up SAEGUS
=============================================

After Anaconda, simuPOP and h5py have been installed, the following steps are used to set up and execute SAEGUS:

+ Clone the 'saegus' repository
+ Install 'saegus'
+ Open a Python terminal, preferrably a Jupyter notebook: ``$ jupyter notebook``
+ Follow the [user guide](https://saegus.readthedocs.io/en/latest/index.html), starting with the [creation of a population from raw data](https://saegus.readthedocs.io/en/latest/population_from_raw_data.html)

Clone the SAEGUS repository
---------------------------------------------

Using HTTPS protocol:

```bash
$ git clone https://github.com/maizeatlas/saegus.git
```

Using SSH protocol

```bash
$ git clone git@github.com:maizeatlas/saegus.git
```

Alternatively if you already have the repository on your machine then use `git pull` command to retrieve the most recent version of the repository. Navigate to the location of the repository on your machine and run:

```bash
$ git pull
```

Install SAEGUS
-----------------------------

In the ``saegus`` directory navigate to the saegus_project/src/ directory. You should see: 

```bash
~/saegus/src$ ls
~/saegus/src$ MANIFEST.in . build . data . docs . saegus . setup.py
```
Run the python setup.py install command to install the saegus package. You should see something like this:

```bash
~/saegus/src$ python setup.py install
running install
running build
running build_py
creating build
creating build/lib
...
```

Open a Jupyter Notebook
-------------------------------------------------------------

In your terminal navigate to the ``$ ~/src/data/`` directory. This is where all the example data files are located for the user guide. In your terminal run:

```bash
$ jupyter notebook
```
A window in your browser will open. This is a much nicer environment to work in than a plain Python interpreter. Especially with saegus' long function names.

Creating a Population from Real Data
------------------------------------------------------------------

Open a new tab in your web browser and follow the walk through: [saegus walkthrough](https://saegus.readthedocs.io/en/latest/index.html)

