#!/bin/bash
echo "`date`"
# INITIAL WARNING: ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD). It is required by scripts running in parallel, and extended to all scripts. 
# Splits libraries according to barcodes specified in mapping file
# Options used:
# -L 1000	Maximum sequence length, in nucleotides [default: 1000]
# -l 250	Minimum sequence length, in nucleotides [default: 200]
# -s 25		Min average qual score allowed in read [default: 25]
# -w 50		Enable sliding window test of quality scores. See split_libraries.py man page for further explanantion.
# -H 13		Maximum length of homopolymer run [default: 6] **Default value is too low for fungi in general**
# -t		Calculate sequence lengths after trimming primers and barcodes [default: False] **Verified not working when coupled with -H option [known bug]**
# -b 0		Barcode type or a number representing the length of the barcode [default: golay_12] **My libraries have been demultiplexed by sequencing facility and provided with no barcode. New identifier added with option --added_demultiplex_field**
# -M 3		Maximum number of primer mismatches [default: 0]
#--reverse_primer_mismatches 3	Set number of allowed mismatches for reverse primers (option -z). [default: 0]
#-z truncate_only	Enable removal of the reverse primer and any subsequence sequence from the end of each read. See split_libraries.py man page for further explanantion. **Actually not needed if pre-process requires denoising**
# --added_demultiplex_field=tag	Adds a field to use in the mapping file as an additional demultiplexing option to the barcode. See split_libraries.py man page for further explanantion.
echo "`date` :: Running split_library.py with forward primer"
split_libraries.py -f $PWD/risinnova_AMF_soil_all_libraries.fasta -q $PWD/risinnova_AMF_soil_all_libraries.qual -o $PWD/risinnova_AMF_soil_left -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -L 1000 -l 250 -s 25 -w 50 -H 13 -t -b 0 -M 3 --reverse_primer_mismatches 3 -z truncate_only --added_demultiplex_field=tag
echo "Done." 
echo "`date` :: Running split_library.py with reverse primer"
split_libraries.py -f $PWD/risinnova_AMF_soil_all_libraries.fasta -q $PWD/risinnova_AMF_soil_all_libraries.qual -o $PWD/risinnova_AMF_soil_right -m $PWD/check_id_map_output/lib_7-8_reverse.csv_corrected.txt -L 1000 -l 250 -s 25 -w 50 -H 13 -t -b 0 -M 3 --reverse_primer_mismatches 3 -z truncate_only --added_demultiplex_field=tag -n 100000
echo "Done."
# This script will denoise a flowgram file in .sff.txt format. It takes a while to complete, depends on libraries dimension.
# left and right stand for forward and reverse, respectivelly.
# The .sff.txt file is the output of sffinfo [sffinfo usage: sffinfo file.sff > file.sff.txt]
# -n 40		Number of CPUs [default: 1]
# --titanium	Select Titanium defaults for denoiser, otherwise use FLX defaults [default: False]
echo "`date` :: Running denoise_wrapper.py on left"
denoise_wrapper.py -v -i $PWD/risinnova_AMF_soil_all_libraries.sff.txt -f $PWD/risinnova_AMF_soil_left/seqs.fna -o $PWD/risinnova_AMF_soil_left/denoised -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -n 40 --titanium
echo "Done." 
echo "`date` :: Running denoise_wrapper.py on right"
denoise_wrapper.py -v -i $PWD/risinnova_AMF_soil_all_libraries.sff.txt -f $PWD/risinnova_AMF_soil_right/seqs.fna -o $PWD/risinnova_AMF_soil_right/denoised -m $PWD/check_id_map_output/lib_7-8_reverse.csv_corrected.txt -n 40 --titanium
echo "Done."
# Inflate denoiser results so they can be passed directly to OTU pickers.
# The inflation process writes each centroid sequence n times, where n is the number of reads that cluster to that centroid, and writes each singleton once. Flowgram identifiers are mapped back to post-split_libraries identifiers in this process (i.e., identifiers in fasta fps)
echo "`date` :: Running inflate_denoiser_output.py on left"
inflate_denoiser_output.py -c $PWD/risinnova_AMF_soil_left/denoised/centroids.fasta -s $PWD/risinnova_AMF_soil_left/denoised/singletons.fasta -f risinnova_AMF_soil_left/seqs.fna -d $PWD/risinnova_AMF_soil_left/denoised/denoiser_mapping.txt -o $PWD/risinnova_AMF_soil_left/denoised/risinnova_AMF_soil_left.inflated.fasta
echo "Done." 
echo "`date` :: Running inflate_denoiser_output.py on right"
inflate_denoiser_output.py -c $PWD/risinnova_AMF_soil_right/denoised/centroids.fasta -s $PWD/risinnova_AMF_soil_right/denoised/singletons.fasta -f $PWD/risinnova_AMF_soil_right/seqs.fna -d $PWD/risinnova_AMF_soil_right/denoised/denoiser_mapping.txt -o $PWD/risinnova_AMF_soil_right/denoised/risinnova_AMF_soil_right.inflated.fasta
echo "Done."
# truncate_reverse_primer takes a demultiplexed fasta file, finds a specified reverse primer sequence, and truncates this primer and subsequent sequences following the reverse primer. This step is required because, eventhou the reverse primer has been truncated during split_libraries, denoising/inflate process writes back the original sequence.
# -M 3	Number of mismatches allowed in the reverse primer. [default: 2]
# (-z, omitted, correctly set to default: truncate_only) 
echo "`date` :: Running truncate_reverse_primer.py on left"
truncate_reverse_primer.py -f $PWD/risinnova_AMF_soil_left/denoised/risinnova_AMF_soil_left.inflated.fasta -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_AMF_soil_left/denoised/ -M 3
echo "Done." 
echo "`date` :: Running truncate_reverse_primer.py on right"
truncate_reverse_primer.py -f $PWD/risinnova_AMF_soil_right/denoised/risinnova_AMF_soil_right.inflated.fasta -m $PWD/check_id_map_output/lib_7-8_reverse.csv_corrected.txt -o $PWD/risinnova_AMF_soil_right/denoised/ -M 3
echo "Done."
# Get the reverse complement of all sequences
# applied on right/reverse sequences
echo "`date` :: Running adjust_seq_orientation.py on right"
adjust_seq_orientation.py -i $PWD/risinnova_AMF_soil_right/denoised/risinnova_AMF_soil_right.inflatedsta_rev_primer_truncated.fna -o $PWD/risinnova_AMF_soil_right/denoised/risinnova_AMF_soil_right.inflatedsta_rev_primer_truncated_rev.fna
echo "Done."
# Merges forward and reverse (but now oriented in same direction) sequences 
echo "`date` :: Merging left and right sequence subsets"
cat $PWD/risinnova_AMF_soil_left/denoised/risinnova_AMF_soil_left.inflatedsta_rev_primer_truncated.fna $PWD/risinnova_AMF_soil_right/denoised/risinnova_AMF_soil_right.inflatedsta_rev_primer_truncated_rev.fna > $PWD/risinnova_AMF_soil_denoised_inflated.fasta
echo "Done."
# Identify chimeric sequences in input FASTA file
# The reference file is the collection of Virtual Taxa Type sequences from MaarjAM database
# -m usearch61	Chimera detection method. Choices: blast_fragments or ChimeraSlayer or usearch61. [default:ChimeraSlayer]
# usearch61 performs both de novo (abundance based) chimera and reference based detection.
# Users can play with usearch61 option to reduce chimeric false positives, but default values do the right job
echo "`date` :: Running identify_chimeric_seqs.py"
identify_chimeric_seqs.py -i $PWD/risinnova_AMF_soil_denoised_inflated.fasta -m usearch61 -o $PWD/usearch61_checked_chimeras/ -r ~/storage/QIIME_DATABASES/maarjAM/vt_types_fasta_from_31-03-2013.txt 
echo "Done."
# Removes sequences from a fasta or fastq file based on input criteria
# Discard all sequences that show up in chimera checking output. NOTE: It is very important to pass -n here as this tells the script to negate the request, or discard all sequences that are listed via -s. This is necessary to remove the identified chimeras from inseqs.fasta.
echo "`date` :: Running filter_fasta.py"
filter_fasta.py -f $PWD/risinnova_AMF_soil_denoised_inflated.fasta -o $PWD/risinnova_AMF_soil_denoised_inflated_filtered.fasta -s $PWD/usearch61_checked_chimeras/chimeras.txt -n
echo "Done."
# An open-reference workflow for OTU picking, taxonomy assignment, phylogenetic tree construction, and OTU table construction, using MaarjAM as the reference collection.
# The open-reference OTU picking process: reads are clustered against a reference sequence collection and any reads which do not hit the reference sequence collection are subsequently clustered de novo.
# taxonomy assignment is performed with blast and MaarjAM database
# phylogenetic tree construction is done by raxml_v730
# Sequence similarity threshold is set to 97% [default: 0.97]
# The alignment used to refine OTU table and to build phylogenetic tree is performed by pyNAST. PyNAST needs a set of aligned reference sequences: vt_types_mafft_aligned.fasta is the alignment of VT Types sequence from MaarjAM, obtained with mafft, with no editing.
#
# Options used:
# -r The reference sequences. The complete MaarjAM database.
# -a Run in parallel where available [default: False]
# -O 40	Number of jobs to start
# -m usearch61	The OTU picking method to use for reference and de novo steps. Passing usearch61 means that usearch61 will be used for the de novo steps and usearch61_ref will be used for reference steps. [default: uclust]
# --min_otu_size 10	The minimum otu size (in number of sequences) to retain the otu [default: 2]
#
# The script behavior is controlled by parameter file.
# 
# OTU picker parameters specified in parameter file
# pick_otus:similarity	0.97
# Multiple sequence alignment parameters
# align_seqs:template_fp	~/storage/QIIME_DATABASES/maarjAM/vt_types_mafft_aligned.fasta
# align_seqs:alignment_method	pynast	
# align_seqs:pairwise_alignment_method	uclust	
# align_seqs:min_length	150
# align_seqs:min_percent_id	75.0
# Alignment filtering (prior to tree-building) parameters
# filter_alignment:allowed_gap_frac	 0.999999
# filter_alignment:remove_outliers	False
# filter_alignment:threshold	3.0
# filter_alignment:suppress_lane_mask_filter	True
# Taxonomy assignment parameters
# assign_taxonomy:id_to_taxonomy_fp	~/storage/QIIME_DATABASES/maarjAM/maarjAM.id_to_taxonomy.5.txt
# assign_taxonomy:reference_seqs_fp	~/storage/QIIME_DATABASES/maarjAM/maarjAM.5.fasta
# assign_taxonomy:assignment_method	blast
# assign_taxonomy:e_value	0.00001
# Phylogenetic tree building parameters
# make_phylogeny:tree_method	raxml_v730
# make_phylogeny:root_method	tree_method_default
echo "`date` :: Running pick_open_reference_otus.py"
pick_open_reference_otus.py -i $PWD/risinnova_AMF_soil_denoised_inflated_filtered.fasta -r ~/storage/QIIME_DATABASES/maarjAM/maarjAM.5.fasta -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10 -p $PWD/parameters_18S_open_reference_usearch61_97.txt -a -O 40 -m usearch61 --min_otu_size 10
echo "Done."
# Conversion of OTU biom file to tabular format for inspection. Please note that VTX code is attached to the species name, with "_" as separator. The file can be easily imported into a spreadsheet and VT isolated for subsequent analyses.  
echo "`date` :: Converting OTU table from biom to tabular format"
biom convert -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures.tab -b --header-key taxonomy
echo "Done."
# Sort the OTU table according to the original sample order. Required to restore the order and get pretty formated output from subsequent steps
echo "`date` :: Sorting OTU table by sample ID"
sort_otu_table.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -l $PWD/sample_order_list.txt
echo "Done."
#Conversion of sorted OTU biom file to tabular format for inspection.
echo "`date` :: Converting sorted OTU table from biom to tabular format"
biom convert -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.tab -b --header-key taxonomy
echo "Done."
# Create an interactive OTU heatmap from an OTU table. The sorted BIOM file is used.
echo "`date` :: Generating OTU heatmap"
make_otu_heatmap_html.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/OTU_HEATMAP
echo "Done." 
# Perform taxonomy summaries and plots. The sorted BIOM file is used.
echo "`date` :: Running summarize_taxa_through_plots.py"
summarize_taxa_through_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/taxa_summary -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt
echo "Done."
# Compile a summary of the information from the OTU table. The sorted BIOM file is used. It will print a summary of the count information on a per-sample basis to the new file specified by the -o parameter. Stats from this file will be used in beta diversity analysis.
echo "`date` :: Summarizing OTU table"
biom summarize-table -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized.txt
echo "Done."
# Summarize an OTU table by a single column in the mapping file. If you have multiple replica for samples, use this script to collapse the replica. Mapping file should be designed accordingly. 
echo "`date` :: Summarizing OTU table by category"
summarize_otu_by_cat.py -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom  -c 'Management_MaturationStage' -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage.biom
biom convert -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage.tab -b --header-key taxonomy
echo "Done."
# Identify the core microbiome. Mapping file should be designed accordingly.
echo "`date` :: Computing core microbiome on selected compartment"
compute_core_microbiome.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/core_microbiome_soil --mapping_fp $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt --valid_states "Compartment:soil"
echo "Done."
# Perform alpha rarefaction. The steps performed by this script are: Generate rarefied OTU tables; compute alpha diversity metrics for each rarefied OTU table; collate alpha diversity results; and generate alpha rarefaction plots.
# -a Run in parallel where available [default: False]
# -O 40	Number of jobs to start
# --num_steps=25	Number of steps (or rarefied OTU table sizes) to make between min and max counts [default: 10] **The higher is this value, the smoother are the curves**
# --max_rare_depth=17000	The upper limit of rarefaction depths [default: median sequence/sample count]
# --min_rare_depth=200	The lower limit of rarefaction depths [default: 10]
# test various --max_rare_depth to obtain pretty figures. Start with default, and increase untill curves (i.e., observed species) reach plateau and are fully contained in the figure.
#
# The script behavior is controlled by parameter file.
#
# Make rarefaction plots parameters
# make_rarefaction_plots:imagetype                png
# make_rarefaction_plots:resolution               400
# make_rarefaction_plots:background_color         white
# multiple_rarefactions:num-reps                  25
# multiple_rarefactions:step                      200
# multiple_rarefactions_even_depth:num-reps       25
# Alpha diversity parameters
# alpha_diversity:metrics	chao1,observed_species,shannon,PD_whole_tree
echo "`date` :: Running alpha_rarefaction.py" 
alpha_rarefaction.py -a -O 40 -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/alpha_rarefaction_SomeMetrics -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt --num_steps=25 --max_rare_depth=17000 --min_rare_depth=200
echo "Done."
# Compute beta diversity distance matrices and generate PCoA plots. Qiime v.1.8 will generate only 3D PCoA Plots.
# -e value (Depth of coverage for even sampling [default: None]) is extracted from the summary file of the sorted OTU table.
#
# The script behavior is controlled by parameter file.
#
# Beta diversity parameters **ALL METRICS**
# beta_diversity:metrics	abund_jaccard,binary_chisq,binary_chord,binary_euclidean,binary_hamming,binary_jaccard,binary_lennon,binary_ochiai,binary_otu_gain,binary_pearson,binary_sorensen_dice,bray_curtis,canberra,chisq,chord,euclidean,gower,hellinger,kulczynski,manhattan,morisita_horn,pearson,soergel,spearman_approx,specprof,unifrac,unifrac_g,unifrac_g_full_tree,unweighted_unifrac,unweighted_unifrac_full_tree,weighted_normalized_unifrac,weighted_unifrac
echo "`date` :: Getting smaller library size"
MIN=`awk '{if ($0 ~ /^ Min/) {split($2,min,"."); print min[1]}}' $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized.txt`
echo "Min library size is: "$MIN
echo "`date` :: beta_diversity_through_plots.py"
beta_diversity_through_plots.py -a -O 40 -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics" -e $MIN -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt
echo "Done."
# Generate 2D PCoA plots using the principal coordinates file generated by performing beta diversity measures of an OTU table. Select the metric you are interested in.
echo "`date` :: make_2d_plots.py -> unweighted_unifrac"
make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/unweighted_unifrac_pc.txt -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/2d_plots
echo "Done."
echo "`date` :: make_2d_plots.py -> weighted_unifrac"
make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/weighted_unifrac_pc.txt -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/2d_plots
echo "Done."
echo "`date` :: make_2d_plots.py -> weighted_normalized_unifrac"
make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/weighted_normalized_unifrac_pc.txt -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics"/2d_plots
echo "Done."
# Perform jackknifed UPGMA clustering and build jackknifed 2d and 3D PCoA plots.
# -e	(Number of sequences to include in each jackknifed subset [REQUIRED]) is set at 75% of smaller sample size.
echo "`date` :: Compute jackknifed_beta_diversity.py with 75% smaller sample size"
let "MIN75=$MIN*75/100"
jackknifed_beta_diversity.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted.biom -m $PWD/check_id_map_output/lib_7-8_forward.csv_corrected.txt -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN75"_AllMetrics_jackknifed" -e $MIN75 -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt
echo "`date` :: ANALYSIS ON SINGLE SAMPLE DONE. STARTING POOLED SAMPLE ANALYSIS"
# The following steps repeats the analysis on the pooled replica.
# Pooled replica BIOM file has to be fixed. It probably contains columns messed up. Fix the columns order in a spreadsheet using the .tab file, and convert it back to BIOM
biom convert -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage.tab -o  $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom --table-type="otu table" --process-obs-metadata taxonomy

biom summarize-table -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited_summarized.txt

make_otu_heatmap_html.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/OTU_HEATMAP_POOLED_REPLICA
# For the following steps a new mapping file has to be done. Mapping file should be designed accordingly to reflect pooled replica composition. Select one replica, and replace the sample name with the pooled replica identifier.
summarize_taxa_through_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/taxa_summary_pooled_replica -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv

compute_core_microbiome.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/core_microbiome_soil_Pooled_Replica --mapping_fp $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv --valid_states "Compartment:soil"

alpha_rarefaction.py -a -O 40 -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/alpha_rarefaction_SomeMetrics_Pooled_Replica -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt --num_steps=25 

MIN=`awk '{if ($0 ~ /^ Min/) {split($2,min,"."); print min[1]}}' $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited_summarized.txt`
echo "Min library size is: "$MIN
beta_diversity_through_plots.py -a -O 40 -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica" -e $MIN -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt

make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/unweighted_unifrac_pc.txt -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/2d_plots

make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/weighted_unifrac_pc.txt -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/2d_plots

make_2d_plots.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/weighted_normalized_unifrac_pc.txt -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN"_AllMetrics_Pooled_Replica"/2d_plots

let "MIN75=$MIN*75/100"
jackknifed_beta_diversity.py -i $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/otu_table_mc10_w_tax_no_pynast_failures_sorted_summarized_by_Management_MaturationStage_edited.biom -m $PWD/check_id_map_output/risinnova_AMF_18S_mapping_pooled_replica.csv -o $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/"beta_div_even"$MIN75"_AllMetrics_Pooled_Replica_jackknifed" -e $MIN75 -t $PWD/risinnova_18S_open_refererence_usearch61_97_VT_filtered_mc10/rep_set.tre -p $PWD/parameters_18S_open_reference_usearch61_97.txt

echo "`date` :: ALL DONE. Bye"

