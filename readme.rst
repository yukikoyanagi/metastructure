This file is a step-by-step guide to metastructure analysis.

Aim
===
Given the primary and secondary structure of a protein, predict beta-strand pairings to deduce the motif.


Software Requirement
====================

- python 3.7
- numpy 1.16
- alignment 1.1.11
- fatgraph 0.0.1


Method -Binary classification
=============================

Prepare the list of proteins
----------------------------
*If you have created new set of Hbond and summary data, remember to modify ``SUM_DIR`` and ``HQ60_DIR`` in ``asfatgraph.py``*
The algorithm currently does not support barrels or bifurcations in proteins. It is therefore necessary to prepare a list of proteins without these structures. Detection of barrels and bifurcations can be done using ``asfatgraph.py``;
::
   
   python asfatgraph.py [pdb id] -f

for detecting bifurcation, and
::
   
   python asfatgraph.py [pdb id] -b

for detecting barrel. Replace ``[pdb id]`` with the 5-letter PDB ID (including chain id). With GNU parallel this can be parallelised. Suppose ``pdbids.txt`` contains all PDB ID's (one per line). Then
::
   
   parallel python asfatgraph.py {} -f :::: pdbids.txt

runs through all proteins in the file. You can then use, for example, ``grep`` to construct a list of all proteins without bifurcation or barrel. Use ``shuff -n200`` to select 200 validation data, and the rest becomes learning data.

Prepare the learning data
-------------------------
The alpha-gamma segments need to be extracted from the learning data. This can be done using ``gammasegs.py``;
::

   parallel python gammasegs.py {} [Dir with summary files] [Dir with Hbond files] -s0 :::: learning.lst

where ``learning.lst`` is the file containing PDB ID of learning data. Vary ``-s`` option from 0 to 4 and save the results to a file (e.g. ``gammasegs.txt``).

Limit the validation data to 20 strands
---------------------------------------
The algorithm is capable of producing results for proteins up to 20 strands. It can probably process larger proteins, but it depends on the partial structure determined by alignment. To find the number of strands, use ``asfatgraph.py`` with ``-n`` switch:
::

   python asfatgraph.py [pdb id] -n

It can of course be combined with GNU parallel. Create a list of validation proteins with <=20 strands (which we call ``validation_20.lst``).

Make partial pairing matrices using alignment
---------------------------------------------
Use ``makepartialmatrix.py``. We build partial matrices up to the 5th diagonal, i.e. by aligning segments with up to 4 beta strands. This can be changed by ``-u`` switch. Run
::
   
   python makepartialmatrix.py -h

to see other options. Note unless you are running this for a few, small proteins, it should be run on a cluster. See ``makepartialmatrix.sh`` for an example. The output is a series of *.mat.pkl files in output_v\* directory. The reduction of partial pairing matrices depending on the strand counts is done in the next step; this could be changed to speed up the processing.


Construct candidate structures
------------------------------
This is done by ``frompartialmatrix.py``. Again for normal runs it should be run on a cluster. See ``frompartialmatrix.sh`` for an example. It produces a series of .pkl files in the specified output directory containing candidate motifs.


Assessment of results
---------------------
*The script requires gbs_arr_20.pkl, which is a (pickled) numpy array for the topology filter.*
Done by ``assess2.py``. Should be run on a cluster. See ``assess.sh`` for an example. Note ``n_strands`` file should have two columns (space separated), the first PDB ID and the second strand cound. The result is a .pkl file with a list of tuples; (PDB ID, size, accepted/rejected, Recall, Precision).


Method2 -Prediction
===================
Create pairing score matrix
---------------------------
Strand pairing matrix of the format ``5VYRA_bp.json`` in sample data folder. This was created by using ``bp2pmat.py`` script from betapro output file [1].

Make prediction
---------------
Done by ``makepreditiction.py``. Use ``-h`` switch to see the options. For an example of running on a cluster, see ``makeprediction.sh``. 
