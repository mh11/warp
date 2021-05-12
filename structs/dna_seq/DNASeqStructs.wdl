version 1.0

struct SampleAndUnmappedBams {
  String base_file_name
  String? final_gvcf_base_name
  Array[String] flowcell_unmapped_bams
  String sample_name
  String unmapped_bam_suffix
}

struct ReferenceFasta {
  String ref_dict
  String ref_fasta
  String ref_fasta_index
  String ref_alt
  String ref_sa
  String ref_amb
  String ref_bwt
  String ref_ann
  String ref_pac
}

struct DNASeqSingleSampleReferences {
  String contamination_sites_ud
  String contamination_sites_bed
  String contamination_sites_mu
  String calling_interval_list

  ReferenceFasta reference_fasta

  Array[String] known_indels_sites_vcfs
  Array[String] known_indels_sites_indices

  String dbsnp_vcf
  String dbsnp_vcf_index

  String evaluation_interval_list

  String haplotype_database_file
}

struct VariantCallingScatterSettings {
   Int haplotype_scatter_count
   Int break_bands_at_multiples_of
}

struct ExomeGermlineSingleSampleOligos {
  String target_interval_list
  String bait_interval_list
  String bait_set_name
}

struct CrossSpeciesContaminationReferences {
  String filter_bwa_image
  String kmer_file
  String meats_bwa_image
  String meats_fasta
  String meats_fasta_dict
  String meats_taxonomy_file
  String microbe_bwa_image
  String microbe_fasta
  String microbe_fasta_dict
  String microbe_taxonomy_file
  String normalization_file
  String metrics_script_file
  Float score_min_identity
  Int reads_after_downsampling
}

struct PapiSettings {
  Int preemptible_tries
  Int agg_preemptible_tries
}
