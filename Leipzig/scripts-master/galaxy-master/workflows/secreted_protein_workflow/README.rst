This is package is a Galaxy workflow for the identification of candidate
secreted proteins from a given protein FASTA file.

It runs SignalP v3.0 (Bendtsen et al. 2004) and selects only proteins with a
strong predicted signal peptide, and then runs TMHMM v2.0 (Krogh et al. 2001)
on those, and selects only proteins without a predicted trans-membrane helix.
This workflow was used in Kikuchi et al. (2011), and is a simplification of
the candidate effector protocol described in Jones et al. (2009).

See http://www.galaxyproject.org for information about the Galaxy Project.


Availability
============

This workflow is available to download and/or install from the main
Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/peterjc/secreted_protein_workflow

Test releases (which should not normally be used) are on the Test Tool Shed:

http://testtoolshed.g2.bx.psu.edu/view/peterjc/secreted_protein_workflow

Development is being done on github here:

https://github.com/peterjc/pico_galaxy/tree/master/workflows/secreted_protein_workflow


Sample Data
===========

This workflow was developed and run on several nematode species. For example,
try the protein set for *Bursaphelenchus xylophilus* (Kikuchi et al. 2011):

ftp://ftp.sanger.ac.uk/pub/pathogens/Bursaphelenchus/xylophilus/Assembly-v1.2/BUX.v1.2.genedb.protein.fa.gz

You can upload this directly into Galaxy via this URL. Galaxy will handle
removing the gzip compression to give you the FASTA protein file which has
18,074 sequences. The expected result (selecting organism type Eukaryote)
is a FASTA protein file of 2,297 predicted secreted protein sequences.


Citation
========

If you use this workflow directly, or a derivative of it, in work leading
to a scientific publication, please cite:

Cock, P.J.A. and Pritchard, L. (2014). Galaxy as a platform for identifying
candidate pathogen effectors. Chapter 1 in "Plant-Pathogen Interactions:
Methods and Protocols (Second Edition)"; P. Birch, J. Jones, and J.I. Bos, eds.
Methods in Molecular Biology. Humana Press, Springer. ISBN 978-1-62703-985-7.
http://www.springer.com/life+sciences/plant+sciences/book/978-1-62703-985-7

Peter J.A. Cock, Björn A. Grüning, Konrad Paszkiewicz and Leighton Pritchard (2013).
Galaxy tools and workflows for sequence analysis with applications
in molecular plant pathology. PeerJ 1:e167
http://dx.doi.org/10.7717/peerj.167

Bendtsen, J.D., Nielsen, H., von Heijne, G., Brunak, S. (2004)
Improved prediction of signal peptides: SignalP 3.0. J Mol Biol 340: 783–95.
http://dx.doi.org/10.1016/j.jmb.2004.05.028

Krogh, A., Larsson, B., von Heijne, G., Sonnhammer, E. (2001)
Predicting transmembrane protein topology with a hidden Markov model:
application to complete genomes. J Mol Biol 305: 567- 580.
http://dx.doi.org/10.1006/jmbi.2000.4315


Additional References
=====================

Kikuchi, T., Cotton, J.A., Dalzell, J.J., Hasegawa. K., et al. (2011)
Genomic insights into the origin of parasitism in the emerging plant
pathogen *Bursaphelenchus xylophilus*. PLoS Pathog 7: e1002219.
http://dx.doi.org/10.1371/journal.ppat.1002219

Jones, J.T., Kumar, A., Pylypenko, L.A., Thirugnanasambandam, A., et al. (2009)
Identification and functional characterization of effectors in expressed
sequence tags from various life cycle stages of the potato cyst nematode
*Globodera pallida*. Mol Plant Pathol 10: 815–28.
http://dx.doi.org/10.1111/j.1364-3703.2009.00585.x


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/peterjc/tmhmm_and_signalp
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id

However, at the time of writing those Galaxy tools have their own
dependencies required for this workflow which require manual
installation (SignalP v3.0 and TMHMM v2.0).


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial release to Tool Shed (May, 2013)
        - Expanded README file to include example data
v0.0.2  - Updated versions of the tools used, inclulding core Galaxy Filter
          tool to avoid warning about new ``header_lines`` parameter.
        - Added link to Tool Shed in the workflow annotation explaining there
          is a README file with sample data, and a requested citation.
======= ======================================================================


Developers
==========

This workflow is under source code control here:

https://github.com/peterjc/pico_galaxy/tree/master/workflows/secreted_protein_workflow

To prepare the tar-ball for uploading to the Tool Shed, I use this:

    $ tar -cf secreted_protein_workflow.tar.gz README.rst repository_dependencies.xml secreted_protein_workflow.ga

Check this,

    $ tar -tzf secreted_protein_workflow.tar.gz 
    README.rst
    repository_dependencies.xml
    secreted_protein_workflow.ga
