**************************
``megalodon_extras merge``
**************************

The ``megalodon_extras merge`` command group contains commands to merge multiple per-read modified base or sequence variant databases.

These commands can assist in deploying Megalodon on an array of compute resources.

-----------------------------------------
``megalodon_extras merge modified_bases``
-----------------------------------------

Merge multiple per-read modified base databases together.

This command contains multi-process capabilities, but may encounter disk I/O bottlenecks.
Note that the full set of modified base positions must be stored in memory to allow this command to process at high performance.
Thus the number of processes should likely be set dependent upon the amount of RAM available and not the number of CPU cores available.
It is recommended that the output location be on a fast disk (i.e. local SSD and not NFS or mounted drives).

-----------------------------------
``megalodon_extras merge variants``
-----------------------------------

Merge multiple per-read sequence variant databases together.
