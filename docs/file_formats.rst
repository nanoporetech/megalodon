************
File Formats
************

This page describes the output file formats produced by ``megalodon``.

------------
Base Calling
------------

Basecalling produces only FASTA format output at this time.

-------
Mapping
-------

Mapped reads can be output in SAM, BAM or CRAM formats.

-----------------------
Per-read Modified Bases
-----------------------

~~~~~~~~
Database
~~~~~~~~

The primary output for per-read modified base results is an `sqlite database <https://www.sqlite.org/index.html>`_.
This database contains an indexed table with per-read, per-position, modified base scores, as well as auxiliary tables with read, modification type and reference position information.
The read table (``read``) contains the read UUID.
The modification type table (``mod``) contains the single letter modified base code, the location sequence match motif, the raw (including ambiguous bases) motif, and the relative modified base position within the motifs.
The reference position table (``pos``) contains the mapped 0-based position, strand (1=forward, -1=reverse) and chromosome (via a final ``chrm`` table which contains the chromosome text).

This database may be accessed via the ``megalodon.mods.ModsDb`` object.

~~~~~~~~~~~~~
Tab-delimited
~~~~~~~~~~~~~

Modified bases results are also available via tab-dellimited text output.
This output can be requested via the ``--write-mods-text`` flag.
This output contains the following fields: ``read_id``, ``chrm``, ``strand``, ``pos``, ``mod_log_prob``, ``can_log_prob``, ``mod_base``, ``motif``

-------------------------
Aggregated Modified Bases
-------------------------

The default aggregated modified base output is the bedMethyl format (`description here <https://www.encodeproject.org/data-standards/wgbs/>`_).
Alternative formats are `wiggle <https://genome.ucsc.edu/goldenPath/help/wiggle.html>`_ (variableStep) and VCF (treating the modified base as if it were a sequence variant).

--------------------------
Per-read Sequence Variants
--------------------------

Docs coming shortly.

----------------------------
Aggregated Sequence Variants
----------------------------

Docs coming shortly.
