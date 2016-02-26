
##The definition of a high-confidence orthologue gene set from separate *de novo* RNA-seq assemblies

This is a step-by-step guide to the bioinformatics analyses that generated a high-confidence orthologue gene set and input for coalescence analyses in the following study:

Nürnberger, Lohse, Fijarczyk, Szymura and Blaxter. Para-allopatry in hybridising fire-bellied toads (*Bombina bombina* and *B. variegata*): inference from transcriptome-wide coalescence analyses. [bioRxiv manuscript](http://biorxiv.org/content/early/2015/10/28/030056)

This repository is provided to document the analyses. The pipeline is highly customised and would need extensive modification in order to be applied to other datasets. At present, there are no plans to turn this pipeline into a software tool. The rationale for the analyses is presented in the manuscript. 


####Input data:
Trinity<sup>1</sup> assemblies of Illumina RNAseq reads (fasta file)

Assemblies (Newbler/Mira<sup>2</sup>/CAP3<sup>3</sup>) of Roche 454 RNA-seq reads (fasta file)

####Additional Dependencies
BLAST, Bowtie<sup>4</sup>, Clustal-omega<sup>5</sup>, MySQL, Muscle<sup>6</sup>,  RSEM<sup>7</sup>, Custal-omega<sup>8</sup>, PAL2NAL<sup>9</sup>, EMBOSS revseq and water, Gblocks<sup>10</sup>


###Import of assembly data into a MySQL database

In the following, one of the four taxa (*B. bombina*) is chosen as an example

A MySQL database is used to store for each taxon information about contigs, Blastx hits to the *Silurana (Xenopus) tropicalis* proteome, predicted open reading frames (ORFs) and read coverage. The latter two analyses are carried out with packages that are bundled with Trinity (see for further information). Here are the commands:

	blastx -db /.../xtropicalis/xtropProtein -query Trinity.fasta -evalue 0.01 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore' -out Bb_blastx_xt.out -num_threads 6 -max_target_seqs 100
	/.../trinity/trinity-plugins/transdecoder/TransDecoder -t Trinity.fasta &
	/.../trinity/trinity-2013-11-10/util/RSEM\_util/run\_RSEM\_align\_n\_estimate --transcripts <fn> --seqType fq --left <(gunzip -c <fn>)  --right <(gunzip -c <fn>) --output_dir <dirname>

MySQL tables: bom\_nmt_contigs, bom\_blastx\_xtrop, bom\_orfs, bom\_trinity\_rsem

###Select the components and contigs that will be analysed in the Paralogue Filter

...i.e. all components that have at least one ORF and more than one read-supported contig.

Create a view that lists all contigs with read support:

	create view zv_bom_trinity_w_rs as select seq_id, seq_name, comp_id from bom_nmt_contigs inner join bom_trinity_rsem using (seq_id) where fpkm > 0;

Determine which components have orfs:

	create view zv_bom_components_with_orf as select comp_id, count(*) as contigs from bom_nmt_contigs inner join bom_orfs using (seq_id) group by comp_id;

and update table bom_components accordingly

	update bom_components, zv_bom_components_with_orf set bom_components.with_orf = zv_bom_components_with_orf.contigs where bom_components.comp_id = zv_bom_components_with_orf.comp_id;

Change all the NULL entries for the with_orf field to 0;

Create a view that identifies all components that have at least one orf and more than one read-supported contig:

	create view zv_bom_multiseq_components as select comp_id, count(*) as contigs from zv_bom_trinity_w_rs inner join bom_components using (comp_id) where with_orf > 0 group by comp_id having contigs > 1;

Also define the refseqs per taxon which are the read-supported contigs with the longest total orf span per component. In previous versions of Trinity, this used to be the contig with seq no = 1, but now this is true in only about 16000/18000 contigs with orfs in the B. bombina assembly. So, this contig needs to be identified and passed on to the crunch analysis as the reference contig.

First, compute  the total orf length per contig

	create view zv_bom_total_orf_length_per_contig as select seq_id, sum(bom_orfs.length) as total_orf_length from zv_bom_trinity_w_rs inner join bom_orfs using (seq_id) group by seq_id;

Second, determine the maximum orf length from the previous query per component

	create view zv_bom_max_orf_length_per_comp as select comp_id, max(total_orf_length) as max_sum_orf from zv_bom_trinity_w_rs inner join zv_bom_total_orf_length_per_contig using (seq_id) group by comp_id;

Finally, determine the identity of the contig that has the longest orf per component:

	create view zv_bom_refseqs as select ct.seq_id, seq_name from zv_bom_trinity_w_rs as ct, bom_orfs as orf, zv_bom_max_orf_length_per_comp as maxorf where ct.seq_id = orf.seq_id and ct.comp_id = maxorf.comp_id and orf.total_orf_length = maxorf.max_sum_orf;

It might be a good idea to write this query to a temporary table ( bom_refseqs ) to speed up the next step.

NB: This query will produce more than one 'reference' contig per component whenever there is more than one contig with the max total orf length per comp. This is true for several hundred components. `Trinity_analyse_comp_batch.pl` can handle this, because it gives reference status to the first 'reference' per comp that it encounters and compares all other contigs to this one. 

Then export all read-supported contigs that belong to the multiseq component defined above and highlight the reference contigs among them.

	select ct.seq_name, ref.seq_id as 'ref' from zv_bom_trinity_w_rs as ct left join zv_bom_refseqs as ref using (seq_id) where comp_id in (select comp_id from zv_bom_multiseq_components) into outfile 'bom_crunch_contigs.txt';

Copy this file to the home directory and run it through `Trinity_crunch_input_file.pl` to convert to old style contig names (compxxx_cxx_seqxx), sort the contigs numerically and recode the seq_ids (field 1) into either 1 (ref) or 0 (other). This is the input file for the crunch analysis.


