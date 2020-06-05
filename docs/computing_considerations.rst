************************
Computing Considerations
************************

This page aims to describe the megalodon processing workflow, highlighting relevant computing considerations.

------------------
Raw Signal Loading
------------------

Raw signal is loaded from either single- or multi-FAST5 format via the ``ont_fast5_api``.
Raw signal is loaded within worker processes when a read_id/file is pulled off of the reads queue.

The reads queue is filled from a separate worker process.
This process will enumerate FAST5 files and all read_ids stored within these files, but per-read processing will begin as soon as the first read_id/file is found.
Users may notice a period of time where the progress bar does not have a known total number of reads.
Once read enumeration is complete the progress bar will update to include the known total and ETA for run completion.

------------
Base Calling
------------

Basecalling is performed by the pyguppy backend.
Basecalling consists of running the neural network and then decoding this output.
See `guppy documentation on the community page (login required) <https://community.nanoporetech.com/protocols/Guppy-protocol>`_ for more details.

-----------------
Reference Mapping
-----------------

Read mapping is completed using the ``minimap2`` python interface (``mappy``).
The reference index is loaded into shared memory.
A separate thread is linked to each per-read processing worker in order to access the shared memory index.
Thus users may notice threads opened for this processing.
These threads will generally consume very little compute.

---------------------------------
Variant and Modified Base Calling
---------------------------------

Sequence variant and modified base calling is computed within the per-read processing workers using CPU resources.
Generally, this portion of processing will comsume a minority of the compute resources.
Proposing many variants (e.g. all possible 2+ base indels) or modified bases in all contexts may show a bottle neck at this portion of processing.
Internal testing shows that proposal of all possible single base substitutions shows minimal processing at this portion of per-read processing.

---------------
Writing to Disk
---------------

As of version 2.0, the status of output queues is displayed by default.
If any of these status bars indicate a full queue, megalodon will stall waiting on that process to write data to disk.
Moving the  ``--output-directory`` to a location with faster disk I/O performance should imporove performance.
