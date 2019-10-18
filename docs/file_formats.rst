************
File Formats
************

This page describes the output file formats produced by ``megalodon``.

------------
Base Calling
------------

Basecalling produces only FASTA format output at this time.
Basecalls will be output into the ``basecalls.fasta`` file within the ``--output-directory``.

Basecall anchored modified base calls are output into a custom HDF5 format (similar to the guppy output format descibed on the community page).
The HDF5 format contains a single dataset, ``mod_long_names`` at the root level which contains the modified base long names as described in the model used for calling.
The ``Reads`` group contains all of the per-read modified base scores.
Within this group each reads modified base scores are stored in a dataset indexed by the read_id.
This dataset contains an array with dimensions ``basecall length`` by ``number of modified bases``.
Unlike the guppy output, scores are only recorded where applicable given the basecall made (e.g. 5mC calls are only output at canonical C basecall positions).
Invalid positions are represented with a ``NAN`` value.
Note that the modified base scores outside of the dependent canonical contexts are not effected during training, so these values should not be used from the guppy output.

-------
Mapping
-------

Mapped reads can be output in SAM, BAM or CRAM formats.
Basecalls will be output into the ``mappings.sam``, ``mappings.bam``, or ``mappings.cram`` file within the ``--output-directory``.

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
