.. module:: ase.calculators.exciting

========
exciting
========

.. image:: ../../static/exciting.png

Introduction
============

``exciting`` is a full-potential *all-electron*
density-functional-theory (DFT) package based on the
linearized augmented planewave (LAPW) method. It can be
applied to all kinds of materials, irrespective of the atomic species
involved, and also allows for the investigation of the core
region. The website is http://exciting-code.org/

The module depends on lxml  https://lxml.de/


There are two ways to construct the exciting calculator.

1. Using keyword arguments to specify groundstate attributes
2. Use paramdict to specify an input file  structure with a data structure of dictionaries.

See also the tutorial on the web page of exciting: http://exciting-code.org/lithium-atomic-simulation-environment


Constructor with Groundstate Keywords
-------------------------------------

One is by giving parameters of the ground state in the
constructor. The possible attributes can be found at
http://exciting-code.org/ref:groundstate::

    Exciting(bin='excitingser', kpts=(4, 4, 4), xctype='LDA_PW')


Parameter Dictionary
--------------------

When the paramdict keyword is used, the calculator translates the dictionary given into the exciting XML file format.
Note $EXCITINGROOT environmental variable should be set: details at http://exciting-code.org/tutorials-boron

.. literalinclude:: exciting.py

The calculator constructure above is used to create this exciting input file:

.. highlight:: xml

::

    <?xml version='1.0' encoding='UTF-8'?>
    <input xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.exciting-code.org/excitinginput.xsd">
      <title>N3O</title>
      <structure speciespath="$EXCITINGROOT/species/" autormt="false">
        <crystal>
          <basevect>1.88972595820018 0.00000000000000 0.00000000000000</basevect>
          <basevect>0.00000000000000 1.88972595820018 0.00000000000000</basevect>
          <basevect>0.00000000000000 0.00000000000000 1.88972595820018</basevect>
        </crystal>
        <species chemicalSymbol="N" speciesfile="N.xml">
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
        </species>
        <species chemicalSymbol="O" speciesfile="O.xml">
          <atom coord="0.50000000000000 0.50000000000000 0.50000000000000"/>
        </species>
      </structure>
      <relax/>
      <groundstate tforce="true" ngridk="1 2 3"/>
      <properties>
        <dos/>
        <bandstructure>
          <plot1d>
            <path steps="100">
              <point coord="0.75000   0.50000   0.25000" label="W"/>
              <point coord="0.50000   0.50000   0.50000" label="L"/>
              <point coord="0.00000   0.00000   0.00000" label="G"/>
              <point coord="0.50000   0.50000   0.00000" label="X"/>
              <point coord="0.75000   0.50000   0.25000" label="W"/>
              <point coord="0.75000   0.37500   0.37500" label="K"/>
            </path>
          </plot1d>
        </bandstructure>
      </properties>
    </input>

The translation follows the following rules:
String values are translated to attributes. Nested dictionaries are translated to sub elements.
A list of dictionaries is translated to a list of sub elements named after the key of which the list is the value.
The special key "text()" results in text content of the enclosing tag.


Muffin Tin Radius
=================

Sometimes it is necessary to specify a fixed muffin tin radius different from the default. The muffin tin radii can be set by adding a custom array to the atoms object with the name "rmt":


.. highlight:: python

::

    atoms.new_array('rmt', np.array([-1.0, -1.0, 2.3, 2.0] * Bohr))


Each entry corresponds to one atom. If the rmt value is negative, the default value is used. This array is correctly updated if the atoms are added or removed.

Exciting Calculator Class
=========================

.. autoclass:: ase.calculators.exciting.Exciting




