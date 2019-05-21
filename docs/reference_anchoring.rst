********************
Reference Aanchoring
********************

Megalodon's functionality centers on the anchoring of high-information neural network basecalling output to a reference sequence. Given anchored neural network output, alternatives to the reference (either modified bases or canonical bases) are proposed and scored to produce the highest accuracy results.

The neural network output is anchored to the reference via standard read mapping of produced basecalls to the reference sequence. If no reference mapping is produced (using ``minimap2`` via the ``mappy`` python interface) that read is not processed further (basecalls will be output if requested). This standard read mapping is processed to produce a matching of each basecall with a reference position. Positions within an insertion or deletion are left justified to the previous valid mapping position. This constitutes the reference anchoring used for modified base and SNP calling steps.
