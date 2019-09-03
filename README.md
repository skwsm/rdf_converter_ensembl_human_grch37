# RDF converter for Ensembl Human GRCh37



## Usage

    ruby rdf_converter_ensembl_grch37.rb -g ensembl75_grch37_01.txt -e mart_export_ensembl75_grch37_exon_01.txt -x mart_export_ensembl75_grch37_ext_01.txt > ensembl75_grch37_01.ttl

### Input files

First argument is a tab-delimited file exported from Ensembl Mart service, whose header line is blow.

    Gene stable ID  Transcript stable ID    Protein stable ID       Exon stable ID  Gene name       Gene description        Transcript name Chromosome/scaffold name        Gene start (bp) Gene end (bp)   Strand  Transcript start (bp)   Transcript end (bp)     Gene type


Second argument is a tab-delimited file exported from Ensembl Mart service, whose header line is blow.

    Gene stable ID  Transcript stable ID    Exon stable ID  Exon region start (bp)  Exon region end (bp)    Exon rank in transcript

Third argument is a tab-delimited file exported from Ensembl Mart service, whose header line is blow.

    Gene stable ID  Transcript stable ID    HGNC ID UniProtKB/Swiss-Prot ID UniProtKB/TrEMBL ID

Ensembl BioMart service for GRCh37 is availele at [here](https://grch37.ensembl.org/biomart/martview).
