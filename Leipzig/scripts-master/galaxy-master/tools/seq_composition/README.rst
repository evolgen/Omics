Galaxy tool reporting sequence composition
==========================================

This tool is copyright 2014 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using Biopython library functions) to
loop over given sequence files (in a range of formats including FASTA, FASTQ,
and SFF), and report the count of each letter (i.e. amino acids or bases).

This can be useful for sanity checking assemblies (e.g. proportion of N
bases) or looking at differences in base composition.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_composition


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``seq_composition.py`` (the Python script)
* ``seq_composition.xml`` (the Galaxy tool definition)

The suggested location is in a dedicated ``tools/seq_composition`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="seq_composition/seq_composition.xml" />

You will also need to install Biopython 1.62 or later. If you want to run the unit
tests, include this line in ``tools_conf.xml.sample`` and the sample FASTA files
under the ``test-data`` directory. Then::

    ./run_functional_tests.sh -id seq_composition

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
======= ======================================================================


Developers
==========

This script and related tools are being developed on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/seq_composition

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf seq_composition.tar.gz tools/seq_composition/README.rst tools/seq_composition/seq_composition.py tools/seq_composition/seq_composition.xml tools/seq_composition/tool_dependencies.xml test-data/four_human_proteins.fasta test-data/four_human_proteins.seq_composition.tabular test-data/ecoli.fastq test-data/ecoli.seq_composition.tabular test-data/MID4_GLZRM4E04_rnd30_frclip.sff test-data/MID4_GLZRM4E04_rnd30_frclip.seq_composition.tabular


Check this worked::

    $ tar -tzf seq_composition.tar.gz
    tools/seq_composition/README.rst
    tools/seq_composition/seq_composition.py
    tools/seq_composition/seq_composition.xml
    tools/seq_composition/tool_dependencies.xml
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.seq_composition.tabular
    test-data/ecoli.fastq
    test-data/ecoli.seq_composition.tabular
    test-data/MID4_GLZRM4E04_rnd30_frclip.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.seq_composition.tabular


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
