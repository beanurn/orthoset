
##The definition of a high-confidence orthologue gene set from separate *de novo* RNA-seq assemblies

This is a step-by-step guide to the bioinformatics analyses that generated a high-confidence orthologue gene set and input for coalescence analyses in the following study:

Nürnberger, Lohse, Fijarczyk, Szymura and Blaxter. Para-allopatry in hybridising fire-bellied toads (*Bombina bombina* and *B. variegata*): inference from transcriptome-wide coalescence analyses. [bioRxiv manuscript](http://biorxiv.org/content/early/2015/10/28/030056)

This repository is provided to document the analyses. The pipeline is highly customised and would need extensive modification before it could be applied to other datasets. The rationale for the analyses is presented in the manuscript. We also provide a MySQL dump of the database that was used for these analyses (`bombina.sql.gz`).


####Input data:
Trinity<sup>1</sup> (v.2014-04-13p1) assemblies of Illumina RNA-seq reads (fasta file)

Assemblies (Newbler/Mira<sup>2</sup>(v.3)/CAP3<sup>3</sup>) of Roche 454 RNA-seq reads (fasta file)

####Dependencies
BLAST (v.2.2.26), Bowtie1<sup>4</sup> (v.1.0.0), EMBOSS<sup>5</sup> revseq, MCL<sup>6</sup> (v14-137), Muscle<sup>7</sup> (v3.8.31), MySQL (v5.5.47), OrthoMCL<sup>8</sup> (v2.0.9), RSEM<sup>9</sup> (v1.2.25).


###Import of assembly data into a MySQL database

Up to the point of the OrthoMCL analysis, the pipeline is presented for one of the taxa, *Bombina bombina*.

A MySQL database is used to store for each taxon information about contigs, Blastx hits to the *Silurana (Xenopus) tropicalis* proteome, predicted open reading frames (ORFs) and read coverage. The latter two analyses are carried out with packages that are bundled with Trinity (see [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki) for further information). Here are the commands:

	blastx -db /.../xtropicalis/xtropProtein -query Trinity.fasta -evalue 0.01 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore' -out Bb_blastx_xt.out -num_threads 6 -max_target_seqs 100
	/.../trinity/trinity-plugins/transdecoder/TransDecoder -t Trinity.fasta &
	/.../trinity/util/align_and_estimate_abundance.pl --transcripts Bbom_trinity_140915_clean.fasta --seqType fq --left <(zcat ../Bbom_cut_trim_paired_1.fq.gz) --right <(zcat ../Bbom_cut_trim_paired_2.fq.gz) --est_method RSEM --aln_method bowtie --max_ins_size 500 --output_dir RSEM_paired --trinity_mode --prep_reference &

MySQL tables: bom\_nmt\_contigs, bom_blastx\_xtrop, bom\_orfs, bom\_trinity\_rsem

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

Third, determine the identity of the contig that has the longest orf per component:

	create view zv_bom_refseqs as select ct.seq_id, seq_name from zv_bom_trinity_w_rs as ct, bom_orfs as orf, zv_bom_max_orf_length_per_comp as maxorf where ct.seq_id = orf.seq_id and ct.comp_id = maxorf.comp_id and orf.total_orf_length = maxorf.max_sum_orf;

It might be a good idea to write this query to a temporary table ( bom_refseqs ) to speed up the next step.

NB: This query will produce more than one 'reference' contig per component whenever there is more than one contig with the max. total ORF length in a component. This is true for several hundred components per assembly. `Trinity_analyse_comp_batch.pl` can handle this, because it gives reference status to the first maxORF contigs that it encounters and compares all other contigs to this one. 

Then export all read-supported contigs that belong to the multiseq component defined above and highlight the reference contigs among them.

	select ct.seq_name, ref.seq_id as 'ref' from zv_bom_trinity_w_rs as ct left join zv_bom_refseqs as ref using (seq_id) where comp_id in (select comp_id from zv_bom_multiseq_components) into outfile 'bom_crunch_contigs.txt';

Copy this file to the home directory and run it through `Trinity_crunch_input_file.pl` to convert to old style contig names (compxxx_cxx_seqxx), sort the contigs numerically and recode the seq_ids (field 1) into either 1 (ref) or 0 (other). This is the input file for the crunch analysis (`bom_crunch_contigs.txt.sorted`).

###Identify paralogue pairs within Trinity components

Place perl files into an analysis directory and create a subdirectory within it called `input_files`. Place into this subdirectory the following files for each taxon (here: B. bombina):
`bom_crunch_contigs.txt.sorted`, `bom_transdecoder.cds` (= TransDecoder output file), `bom_blastx_xt.out`.
Also, run `RSEM_results_sorted_by_comp_numerical.pl` to get a numerically sorted version of `RSEM.isoforms.results` (`bom_RSEM.isoforms.results.sorted`) and place this into `input_files` as well.

For the pairwise analysis of contigs within components run `Trinity_analyse_comp_batch.pl`. Make sure the path to MUSCLE is specified. For each pair, `Trinity_analyse_comp_batch.pl` calls `Trinity_align_components4.pl` which does a MUSCLE alignment of the contigs, defragments the alignment as required by calling `Trinity_defragment_alignments_2seq.pl` and collects BLASTX and ORFs annotations from the input files. These are appended to the alignment file (.afa). The names of these alignments are written to `aligned_sequences.txt`. Every 50 alignments, this file is passed to `Trinity_analyse_components2.pl` which scans the alignments for biological feature and computes the sequence similarity in good and bad alignment stretches. These data are appended to `bom_component_crunch.txt`. 

The command line for `Trinity_analyse_comp_batch.pl`is 

	perl Trinity_analyse_comp_batch.pl <fn> <start_comp> <stop_comp>

where fn is `bom_crunch_contigs.txt.sorted`. The same taxon identifier ('bom') should be used for all associated input files. start_comp and end_comp define the subset of the components to be analysed. Choose the first component as start and a non-existent end component to run the analysis for all components. The output file `bom_component_crunch.txt` contains a header line which is repeated every so often to facilitate on-screen viewing of the file. To prepare these data for import into MySQL run it through `mysqlprep_component_crunch.pl` which produces `bom_component_crunch_ed.txt`. Upload this table into MySQL. The commands for the creation of the associated table and for the data import are a bit unwieldy and can be found here: `crunch_table_mysql_commands.txt`. Paralogue status is recorded in the 'category' field of that table like so:

	UPDATE bom_component_crunch SET  category = 'para' WHERE goodMean < 0.98;
	UPDATE bom_component_crunch SET  category = 'para' WHERE goodMean >= 0.98 and poorMean > 0.5;

Export reference contigs and paralogues as follows. Create a temporary table of reference sequences as above (`bom_refseqs`):

	create temporary table bom_paralogues as select concat (component, '_', seqno) as seq_name from bom_component_crunch where category = 'para' and seq_id not in (select seq_id from bom_refseqs);
	select distinct seq_name from bom_paralogues into outfile 'bom_paralogues.txt';


###Identify paralogues among contigs from the Roche 454 assembly

First generate a BLAST database from the Trinity reference contigs (see previous section). Use BLASTN to search this database with the Roche 454 contigs:

	blastn -db /.../bom_trinity_db -query bom_nm_input.fasta -out bom454_blastn_trinity.out -num_threads 6 -outfmt '6 std qlen slen' -evalue 1e-5 -max_target_seqs 10 &

Create an input file for the new paralogue search from the BLAST output file with `Trinity454_define_new_components.pl`. This writes new "components", each consisting of one Roche 454 contig and its BLAST hits (max. 10), to a file that serves as inout for the next round of pairwise sequence comparisons. It also records the highest scoring match for each Roche 454 contig in a separate file. Use `Trinity454_analyse_comp_batch.pl`, `Trinity454_align_components.pl` and `Trinity454_analyse_components.pl` to do the pairwise comparisons. Contigs, BLASTX hits to *S. tropicalis* and ORF predictions need to be provided in a single file each (containing Trinity and Roche 454 data/annotations). No coverage information is required for this step. Process the output file `component_crunch_bom454.txt`, upload data to the MySQL database and identify paralogues as before. Note that the fields are slightly different than before.

The Roche 454 contigs are partitioned into a paralogue and a non-paralogue set. To export the former use

	select distinct component, seqno, seq_name from bom_component_crunch454 where category = 'para' into outfile 'bom_454_paralogues.txt';  # these are the Roche 454 paralogues and their Trinity 'paralogue partners'

The latter are all Roche 454 contigs that have no paralogue status in any pairwise comparison:

	create temporary table bom_454_paralogues as select distinct seq_name from bom_component_crunch454 where category = 'para' and seqno = 'i0';  # all Roche 454 contigs with paralogue status in at least one comparison
	select distinct seq_name from bom_component_crunch454 where seqno = 'i0' and seq_name not in (select seq_name from bom_454_paralogues) into outfile 'bom_454_refseqs.txt';

Each 454 paralogue is given a new name which 'adds' it to the Trinity component of its paralogue partner (the best BLASTN hit if possible). This is done with `add_454paralogues_to_trinity_components.pl`, which reads in the MySQL output file of paralogue pairs (`bom_454_paralogues.txt`), the file of best BLASTN hits generated by `Trinity454_define_new_components.pl` (see above) and a fasta file of Roche 454 contigs with paralogue status. Sequence numbers for these new component members start with 100 (> than the largest sequence number per component in the Trinity assemblies). The script outputs the fasta file with new contigs names and a translation file `bom454_paralogues_names.txt` which records pairs of new and old names.

###OrthoMCL prep

Per taxon the contigset is a concatenation of the following subsets: Trinity refseqs (from file `bom_crunch_contigs.txt`), Trinity paralogues and, if available, Roche 454 refseqs (not really a good term, these are the 'non-paralogues') and Roche 454 paralogues. For the OrthoMCL analysis, it is very important that there are no duplicate entries for any contig (the programme will crash otherwise). Use `check_no_dups_in_seqset.pl` to make sure. Run the TransDecoder analysis as before to identify ORFs.

For contigs with ORFS in plus and minus orientation, only the orientation with the longer total ORF length was kept, because it would not have been desirable for our analysis to obtain more than one cluster from any one set of orthologous contigs. The 'majority orientation' is recorded in the component_crunch tables in the MySQL database and extracted as follows:

	select distinct CONCAT ('component','_','seqno') as contig, orientation from bom_component_crunch into outfile 'bom_trinity_orientation.txt';
	select distinct select seq_name, 454dir, orientation from bom_component_crunch454 into outfile 'bom_454_orientation';

These two files and `bom454_paralogues_names.txt` (see above) form the input for `collate_orientation_info.pl`, which determines the majority orientation either from the input or computes it de novo if necessary. It also replaces the names of Roche 454 paralogues with their new Trinity-like names which are used in the OrthoMCL analysis. Feed the output file, together with the TransDecoder protein file (`.pep`) and a three letter taxon prefix into `drop_false_orfs.pl` to eliminate all ORFs with the wrong orientation. The contig names in the output file have the taxon prefix and are in a format acceptable by OrthoMCL. Also, the stop codon markers (`*`) are removed from the peptide sequences. This file is now ready for OrthoMCL.


###OrthoMCL analysis

The OrthoMCL analysis was carried out with four *Bombina* taxa and the *S. tropicalis* proteome according to the [User guide](http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt).

All-versus-all BLASTP analysis:

	blastp -db goodProteins.fasta -query goodProteins.fasta -out bombina_blast.out -num_threads 6 -outfmt 6 -evalue 1e-5 -seg yes -max_target_seqs 100000

The final MCL clustering step was run with three inflation values: 0.5, 1.5, 5.0. Here are the commands for 1.5:

	/.../orthomcl/mcl-14-137/bin/mcl mclInput --abc -I 1.5 -o mclOutput_1.5
	/.../orthomcl/orthomcl-2.0.9/bin/orthomclMclToGroups G_ 1 < mclOutput_1.5 > groups_1.5.txt

###OrthoMCL cluster filtering

Use `groups_1.5.txt` as input for `tabulate_orthogroups.pl` which gives a detailed listing of the composition of each cluster (`groups_1.5.tab`) and outputs the subset of clusters that fulfill the following three criteria (`groups_1.5_GoodGroups.txt`):
* at least one contig from each *Bombina* taxon
* no paralogue pair from any *Bombina* taxon
* no more than three *S. tropicalis* contigs

`robust_good_groups.pl` takes the `groups_1.5_GoodGroups.txt` file together with the other two MCL output files (`groups_#.#.txt`) as input and determines the subset of 'good' clusters that were formed identically across all three inflation values: `RobustGroups.txt`. These robust orthogroups are used for pairwise sequence alignments between *Bombina* taxa. 


###References

1. Grabherr, M. G., B. J. Haas, M. Yassour, J. Z. Levin, D. A. Thompson, I. Amit, X. Adiconis, L. Fan, R. Raychowdhury, Q. Zeng, Z. Chen, E. Mauceli, N. Hacohen, A. Gnirke, N. Rhind, F. di Palma, B. W. Birren, C. Nusbaum, K. Lindblad-Toh, N. Friedman, and A. Regev. 2011. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat. Biotechnol. 29:644–652.

2. Chevreux, B., T. Pfisterer, B. Drescher, A. J. Driesel, W. E. G. Müller, T. Wetter, and S. Suhai. 2004. Using the miraEST Assembler for Reliable and Automated mRNA Transcript Assembly and SNP Detection in Sequenced ESTs. Genome Res. 14:1147–1159.

3. Huang, X., and A. Madan. 1999. CAP3: a DNA sequence assembly program. Genome Res. 12:868–877.

4. Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10:R25. 

5. Rice, P., I. Longden, and A. Bleasby. 2000. EMBOSS: The European Molecular Biology Open Software Suite. Trends Genet. 16:276–277.

6. Stijn van Dongen, Graph Clustering by Flow Simulation.  PhD thesis, University of Utrecht, May 2000.

7. Edgar, R. C. 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 32:1792–1797.

8. Li, L., C. J. Stoeckert, and D. S. Roos. 2003. OrthoMCL: Identification of Ortholog Groups for Eukaryotic Genomes. Genome Res. 13:2178–2189.

9. Li, B., and C. N. Dewey. 2011. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12:323.










