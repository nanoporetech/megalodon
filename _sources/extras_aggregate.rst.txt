******************************
``megalodon_extras aggregate``
******************************

The ``megalodon_extras aggregate`` command group contains a single command, ``run``, to perform aggregation of per-read sequence variant or modified base results.

----------------------------------
``megalodon_extras aggregate run``
----------------------------------

Aggregate per-read sequence variant and/or modified base from the main ``megalodon`` command.

This command can be useful in processing Megalodon pipelines efficiently.
This command allows the ``megalodon`` command can be performed on one set of computing resources and then ``megalodon_extras aggregate run`` can be completed on a separate set of computing resources.
The ``megalodon`` command, running the basecalling backend, generally requires GPU resources, while the aggregation step generally requires fast disk (SSDs) and a lot of CPU cores.
This command allows one to perform these steps separately and on appropriate compute resources.

Additionally, this command can allow for the adjustment of aggregation parameters without the need to repeat the compute expensive basecalling step.
