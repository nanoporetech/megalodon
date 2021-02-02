************************
Computing Considerations
************************

This page aims to describe the Megalodon processing workflow, highlighting relevant computing considerations.

------------------
Raw Signal Loading
------------------

Raw signal is loaded from either single- or multi-FAST5 format via the ``ont_fast5_api``.
Raw signal is loaded within a single process and distributed out the worker processes.
The input queue status bar indicates how many read signals have been loaded and are awaiting processing.
If this status bar is often empty, raw signal extraction from FAST5 files is likely a processing bottleneck.

The reads queue is filled from a separate worker process.
This process will enumerate FAST5 files and all read_ids stored within these files, but per-read processing will begin as soon as the first read_id/file is found.
Users may notice a period of time where the progress bar does not have a known total number of reads.
Once read enumeration is complete the progress bar will update to include the total number of reads found and ETA for run completion.

------------
Base Calling
------------

Basecalling is performed by the pyguppy backend.
Basecalling consists of running the neural network and then decoding this output.
See `guppy documentation on the community page (login required) <https://community.nanoporetech.com/protocols/Guppy-protocol>`_ for more details.
Parameters can be passed directly to the Guppy server initialization call via the ``--guppy-params`` argument.

-----------------
Reference Mapping
-----------------

Read mapping is completed using the ``minimap2`` python interface (``mappy``).
The reference index is loaded into shared memory.
A separate thread is linked to each per-read processing worker in order to access the shared memory index.
Thus users may notice threads opened for this processing.
These threads will generally consume less compute than the worker processes.

---------------------------------
Variant and Modified Base Calling
---------------------------------

Sequence variant and modified base calling is computed within the per-read processing workers using CPU resources.
Generally, this portion of processing will consume a minority of the compute resources.
Proposing many variants (e.g. all possible 2+ base indels) or modified bases in all contexts may show a bottle neck at this portion of processing.
Internal testing shows that proposal of all possible single base substitutions shows minimal processing at this portion of per-read processing.

---------------
Writing to Disk
---------------

As of version 2.0, the status of output queues is displayed by default.
As of version 2.2, the status of the input signal extraction queue is also displayed.
If any of the output status bars indicate a full queue, Megalodon will stall waiting on that process to write data to disk.
if the input signal extraction quque is often empty, raw signal extraction from FAST5 files is likely a processing bottleneck.
Moving the input data or  ``--output-directory`` respectively to a location with faster disk I/O performance should improve performance.
