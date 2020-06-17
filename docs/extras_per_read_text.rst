**********************************
``megalodon_extras per_read_text``
**********************************

The ``megalodon_extras per_read_text`` command group contains commands to convert per-read modified base database statistics to text files.
These files will be TSV files with headers describing the fields contained within the file.

Note that these scripts are single threaded and can be quite slow for reasonable sized runs.

-------------------------------------------------
``megalodon_extras per_read_text modified_bases``
-------------------------------------------------

Extract text format per-read modified base scores from a Megalodon per-read modified base database.

-------------------------------------------
``megalodon_extras per_read_text variants``
-------------------------------------------

Extract text format per-read sequence variant scores from a Megalodon per-read sequence variant database.
