This package contains Galaxy wrappers for a selection of standalone command
line protein analysis tools:

* SignalP 3.0, THMHMM 2.0, Promoter 2.0 from the Center for Biological
  Sequence Analysis at the Technical University of Denmark,
  http://www.cbs.dtu.dk/cbs/

* WoLF PSORT v0.2 from http://wolfpsort.org/

* PSORTb v3 from http://www.psort.org/downloads/index.html

Also, the RXLR motif tool uses SignalP 3.0 and HMMER 2.3.2 internally.

To use these Galaxy wrappers you must first install the command line tools.
At the time of writing they are all free for academic use, or open source.

These wrappers are copyright 2010-2013 by Peter Cock, James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
Contributions/revisions copyright 2011 Konrad Paszkiewicz. All rights reserved.
See the included LICENCE file for details (MIT open source licence).

The wrappers are available from the Galaxy Tool Shed
http://toolshed.g2.bx.psu.edu/view/peterjc/tmhmm_and_signalp 

Citation
========

If you use any of these Galaxy tools in work leading to a scientific
publication, in addition to citing the invididual underlying tools, please cite:

Peter Cock, Bjoern Gruening, Konrad Paszkiewicz and Leighton Pritchard (2013).
Galaxy tools and workflows for sequence analysis with applications
in molecular plant pathology. PeerJ 1:e167
http://dx.doi.org/10.7717/peerj.167

Full reference information is included in the help text for each tool.


Requirements
============

First install those command line tools you wish to use the wrappers for:

1. Install the command line version of SignalP 3.0 and ensure ``signalp`` is
   on the PATH, see: http://www.cbs.dtu.dk/services/SignalP/

2. Install the command line version of TMHMM 2.0 and ensure ``tmhmm`` is on
   the PATH, see: http://www.cbs.dtu.dk/services/TMHMM/

3. Install the command line version of Promoter 2.0 and ensure ``promoter`` is
   on the PATH, see: http://www.cbs.dtu.dk/services/Promoter

4. Install the WoLF PSORT v0.2 package, and ensure ``runWolfPsortSummary``
   is on the PATH (we use an extra wrapper script to change to the WoLF PSORT
   directory, run runWolfPsortSummary, and then change back to the original
   directory), see: http://wolfpsort.org/WoLFPSORT_package/version0.2/

5. Install hmmsearch from HMMER 2.3.2 (the last stable release of HMMER 2)
   but put it on the path under the name ``hmmsearch2`` (allowing it to
   co-exist with HMMER 3), or edit ``rlxr_motif.py`` accordingly.

Verify each of the tools is installed and working from the command line
(when logged in as the Galaxy user if appropriate).


Manual Installation
===================

1. Create a folder ``tools/protein_analysis`` under your Galaxy installation.
   This folder name is not critical, and can be changed if desired - you
   must update the paths used in ``tool_conf.xml`` to match.

2. Copy/move the following files (from this archive) there:

   * tmhmm2.xml (Galaxy tool definition)
   * tmhmm2.py (Python wrapper script)

   * signalp3.xml (Galaxy tool definition)
   * signalp3.py (Python wrapper script)

   * promoter2.xml (Galaxy tool definition)
   * promoter2.py (Python wrapper script)

   * psortb.xml (Galaxy tool definition)
   * psortb.py (Python wrapper script)

   * wolf_psort.xml (Galaxy tool definition)
   * wolf_psort.py (Python wrapper script)

   * rxlr_motifs.xml (Galaxy tool definition)
   * rxlr_motifs.py (Python script)

   * seq_analysis_utils.py (shared Python code)
   * LICENCE
   * README.rst (this file)

3. Edit your Galaxy conjuration file tool_conf.xml (to use the tools) AND
   also ``tool_conf.xml.sample`` (to run the tests) to include the new tools
   by adding::

    <section name="Protein sequence analysis" id="protein_analysis">
      <tool file="protein_analysis/tmhmm2.xml" />
      <tool file="protein_analysis/signalp3.xml" />
      <tool file="protein_analysis/psortb.xml" />
      <tool file="protein_analysis/wolf_psort.xml" />
      <tool file="protein_analysis/rxlr_motifs.xml" />
    </section>
    <section name="Nucleotide sequence analysis" id="nucleotide_analysis">
      <tool file="protein_analysis/promoter2.xml" />
    </section>

   Leave out the lines for any tools you do not wish to use in Galaxy.

4. Copy/move the test-data files (from this archive) to Galaxy's
   subfolder test-data.

5. Run the Galaxy functional tests for these new wrappers with::

    $ ./run_functional_tests.sh -id tmhmm2
    $ ./run_functional_tests.sh -id signalp3
    $ ./run_functional_tests.sh -id Psortb
    $ ./run_functional_tests.sh -id rxlr_motifs

   Alternatively, this should work (assuming you left the name and id as shown
   above in your XML file ``tool_conf.xml.sample``)::

    $ ./run_functional_tests.sh -sid Protein_sequence_analysis-protein_analysis

   To check the section ID expected, use ./run_functional_tests.sh -list

6. Restart Galaxy and check the new tools are shown and work.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial release
v0.0.2  - Corrected some typos in the help text
        - Renamed test output file to use Galaxy convention of *.tabular
v0.0.3  - Check for tmhmm2 silent failures (no output)
        - Additional unit tests
v0.0.4  - Ignore comment lines in tmhmm2 output.
v0.0.5  - Explicitly request tmhmm short output (may not be the default)
v0.0.6  - Improvement to how sub-jobs are run (should be faster)
v0.0.7  - Change SignalP default truncation from 60 to 70 to match the
          SignalP webservice.
v0.0.8  - Added WoLF PSORT wrapper to the suite.
v0.0.9  - Added our RXLR motifs tool to the suite.
v0.1.0  - Added Promoter 2.0 wrapper (similar to SignalP & TMHMM wrappers)
        - Support Galaxy's <parallelism> tag for SignalP, TMHMM & Promoter
v0.1.1  - Fixed an error in the header of the tabular output from Promoter
v0.1.2  - Use the new <stdio> settings in the XML wrappers to catch errors
        - Use SGE style $NSLOTS for thread count (otherwise default to 4)
v0.1.3  - Added missing file whisson_et_al_rxlr_eer_cropped.hmm to Tool Shed
v0.2.0  - Added PSORTb wrapper to the suite, based on earlier work
          contributed by Konrad Paszkiewicz.
v0.2.1  - Use a script to create the Tool Shed tar-ball (removed some stray
          files accidentally included previously via a wildcard).
v0.2.2  - Include missing test files.
v0.2.3  - Added unit tests for WoLF PSORT.
v0.2.4  - Added unit tests for Promoter 2
v0.2.5  - Link to Tool Shed added to help text and this documentation.
        - More unit tests.
        - Fixed bug with RXLR tool and empty FASTA files.
        - Fixed typo in the RXLR tool help text.
        - Updated citation information (Cock et al. 2013).
        - Adopted standard MIT licence.
        - Use reStructuredText for this README file.
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.2.6  - Use the new $GALAXY_SLOTS environment variable for thread count.
        - Updated the suite_config.xml file (overdue).
======= ======================================================================


Developers
==========

This script and other tools were initially developed on the following hg branches:
http://bitbucket.org/peterjc/galaxy-central/src/seq_analysis
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools

For making the "Galaxy Tool Shed" http://community.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ ./tools/protein_analysis/make_tmhmm_and_signalp.sh

This simplifies ensuring a consistent set of files is bundled each time,
including all the relevant test files.
