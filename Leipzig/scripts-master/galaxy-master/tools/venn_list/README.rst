Galaxy tool to draw a Venn Diagram with up to 3 sets
====================================================

This tool is copyright 2011 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using both the Galaxy and Biopython library
functions) to extract ID lists from tabular, FASTA, FASTQ or SFF files to build
sets, which are then drawn using the R limma package function vennDiagram
(called from Python using rpy).

There are just two files to install:

* ``venn_list.py`` (the Python script)
* ``venn_list.xml`` (the Galaxy tool definition)

The suggested location is in the Galaxy folder ``tools/plotting`` next to other
graph drawing tools, or a dedicated ``tools/venn_list`` directory.

You will also need to install Biopython 1.54 or later, and the R/Bioconductor
pacakge limma. You should already have rpy installed for other Galaxy tools.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. The suggested location is in the "Graph/Display Data" section. Simply add
the line::

  <tool file="venn_list/venn_list.xml" />

If you wish to run the unit tests, also add this to ``tools_conf.xml.sample`` and
move/copy the ``test-data`` files under Galaxy's ``test-data`` folder. Then::

    ./run_functional_tests.sh -id venn_list


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.3  - Initial public release.
v0.0.4  - Ignore blank lines when loading IDs from tabular files
v0.0.5  - Explicit Galaxy error handling of return codes
v0.0.6  - Added unit tests.
        - Use reStructuredText for this README file.
        - Adopt standard MIT licence.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.0.7  - Renamed folder and README file.
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf venn_list.tar.gz tools/venn_list/README.rst tools/venn_list/venn_list.* test-data/venn_list1.pdf test-data/venn_list.tabular test-data/rhodopsin_proteins.fasta

Check this worked::

    $ tar -tzf venn_list.tar.gz
    tools/venn_list/README.rst
    tools/venn_list/venn_list.py
    tools/venn_list/venn_list.xml
    test-data/venn_list1.pdf
    test-data/venn_list.tabular
    test-data/rhodopsin_proteins.fasta


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
