Cistrome ChiLin ver. 2.0, db 1.0
=================================

installing
========================
To install the tool, you will first have to cp chilin.conf to chilin.conf.filled
and modify chilin.conf to define site-wide defaults for your system.

In chilin.conf, there are three main sections-- tools, references, and params.

In the tools section, you will have to define the absolute paths to the various
sub-tools or modules in the analysis pipeline.

In the references section, you will have to define the absolute paths to the
reference files (e.g. assembly, DHS file, gene regions, etc.) for the species
you wish to support.  Please read 'Generating Species References' below for
information on how to generate files for your species of interest.

In the params section, you will have the ability to fine-tune and define
site-wide defaults for some of the sub-tools/modules.  NOTE: the defaults
should work for almost all installations--only make modifications only if you
are sure!

After tailoring chilin.conf to suit your system's needs, simple type:
python setup.py install

NOTE: if you are installing system-wide, you may want to add 'sudo' in front
of that command

Generating Species References
=============================

Finally, don't forget to add your species to chilin.conf.


Detailed Document
=============================

cistrome chilin document: http://cistrome.org/chilin
