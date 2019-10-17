This report is the output of a Snakemake pipeline which processes bacterial WGS data. 
Raw sequencing data was trimmed using `Fastp`_ and summarised using `MultiQC`_.
Species detection was performed using `Kraken2`_. 
Assembly was performed using `SKESA`_.
Assembly quality control was performed using `Quast`_ and summarised using `MultiQC`_.
Coverage was assessed mapping reads back to the assembly using `Minimap2`_, followed by processing with `Samtools`_ sort and `Bedtools`_ genomecov, all in one pipe.
MLST was predicted using `mlst`_.
Resistance genes were predicted using the `ABRicate`_ with the `NCBI resistance gene database`_.
Tool versions can be checked at the bottom of this report, in the `Rules`_ section.

.. _Fastp: https://github.com/OpenGene/fastp
.. _MultiQC: https://multiqc.info/
.. _Kraken2: https://ccb.jhu.edu/software/kraken2/
.. _SKESA: https://github.com/ncbi/SKESA
.. _Quast: https://github.com/ablab/quast
.. _Minimap2: https://github.com/lh3/minimap2
.. _Samtools: http://www.htslib.org/
.. _Bedtools: https://bedtools.readthedocs.io/en/latest/
.. _mlst: https://github.com/tseemann/mlst
.. _ABRicate: https://github.com/tseemann/abricate
.. _NCBI resistance gene database: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047

