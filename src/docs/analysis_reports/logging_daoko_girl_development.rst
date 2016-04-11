
.. run-daoko-girl:

======================================
Log of Development for Run: Daoko.girl
======================================

.. daoko-girl-parameters:

Location and Definition of Parameters
=====================================

Genetic Map
~~~~~~~~~~~

A ``genetic_map`` determines how many chromosomes and how many loci are on each
chromosome. The ``genetic_map`` also defines the rates of recombination between
sequential loci as defined by simuPOP_.

.. _simuPOP: http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec5.html#recombination-operator-recombinator


Synthesis Parameters
~~~~~~~~~~~~~~~~~~~~

The purpose of the ``synthesis_parameters`` is to tell the user exactly how the population
was created. The contents of ``synthesis_parameters`` should be able to tell you
where the genotypes of the founders came from and how the founders were mated with each
other. The following items are defined in ``synthesis_parameters``:

   * ``prefounder_file_name``
   * ``prefounder_names``
   * ``founders``
   * ``integer_to_snp``
   * ``snp_to_integer``
   * ``mating_scheme``
   * ``development_population_size``

Trait Parameters
~~~~~~~~~~~~~~~~

``saegus`` was developed with the intention of simulating quantitative phenotypes.
``trait_parameters`` should tell the user exactly how the phenotypes, ``p``, were
calculated. At present ``trait_parameters`` has the following entries:

   * ``number_of_qtl``
   * ``qtl``
   * ``allele_effect_distribution``
   * ``allele_effect_parameters``
   * ``allele_effects``
   * ``heritability``
   * ``epsilon``

Alleles
~~~~~~~

A separate file for the alleles at each locus is being kept in ``alleles_at_1478_loci.hdf``
for the time being.

Analysis Parameters
~~~~~~~~~~~~~~~~~~~

The ability to map the output from TASSEL to the input from ``saegus`` is crucial.
The ``analysis_parameters`` record information which allows the user to make the
connection between input and output.

   * ``population_name``
   * ``scenario``
   * ``generations``
   * ``operating_population_size``
   * ``output_prefix``
   * ``sample_size``
   * ``sample_ind_ids``
   * ``sample_allele_frequencies``
   * ``segregating_loci``
   * ``saegus_to_tassel_loci``
   * ``tassel_to_saegus_loci``

Intermediate Data
~~~~~~~~~~~~~~~~~

There is a great deal of data and parameters moving around inside even a single
simulation of ``saegus``. It is wasteful to generate new data each time I work on
``saegus``. ``intermediate_data`` records parameters or data that in the case
I take a break from developing. ``intermediate_data`` can also be used for debugging
and metadata.

.. code-block:: python

   intermediate_data = shelve.open('daoko_girl_debug_data')
   intermediate_data['allele_frequencies'] = af
   intermediate_data['segregating_allele_frequencies'] = segregating_frame
   intermediate_data['g'] = np.array(magic1478.indInfo('g'))
   intermediate_data['p'] = np.array(magic1478.indInfo('p'))
   intermediate_data['segregating_loci'] = segregating_loci
   intermediate_data['run_name'] = 'daoko_girl'
   intermediate_data.close()



Running the Simulation and Formatting Output for TASSEL
=======================================================

.. code-block:: shell

   $ pwd
   C:\Users\DoubleDanks\wisser\code\rjwlab-scripts\saegus_project\devel\magic\1478

I loaded the **standard** population with 1478 loci made from prefounders 1 through 8.

.. code-block:: python

   magic1478 = sim.loadPopulation('populations\\magic_1478.pop')
   magic1478.dvars()

   {'rep': 0, 'gen': 3}


I loaded the standard genetic map of 1478 loci with recombination rates of 0.01
between sequential loci. The genetic map is a table with the following columns:

   * chromosome
   * relative locus
   * recom_rate

The index of the table serves as the absolute locus.

.. code-block:: python

   genetic_map = pd.read_hdf('parameters\\genetic_map_1478.hdf)
   genetic_map

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>chromosome</th>
          <th>relative_locus</th>
          <th>recom_rate</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>0</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>1</th>
          <td>1</td>
          <td>1</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1</td>
          <td>2</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1</td>
          <td>3</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1</td>
          <td>4</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>5</th>
          <td>1</td>
          <td>5</td>
          <td>0.01</td>
        </tr>
        <tr>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        </tr>
    </table>
    <p>1478 rows × 3 columns</p>
 </div>

