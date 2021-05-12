version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for QC of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  input {
    String input_bam
    String metrics_filename
    Int preemptible_tries
  }

  String flag_file = metrics_filename + ".done"  # finished flag
  Int memory_size = 3 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}
    touch "~{flag_file}"  # flag successful finish
  fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    preemptible: preemptible_tries
  }
  output {
    String quality_yield_metrics = "~{metrics_filename}"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  input {
    String input_bam
    String output_bam_prefix
    Int preemptible_tries
  }

  String flag_file = output_bam_prefix + ".done"  # finished flag
  
  Int memory_size = 7 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=ALL_READS

    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.insert_size_histogram.pdf

    touch "~{flag_file}"  # flag successful finish
  fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    cpu: "1"
  }
  output {
    String base_distribution_by_cycle_pdf = "~{output_bam_prefix}.base_distribution_by_cycle.pdf"
    String base_distribution_by_cycle_metrics = "~{output_bam_prefix}.base_distribution_by_cycle_metrics"
    String insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    String insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    String quality_by_cycle_pdf = "~{output_bam_prefix}.quality_by_cycle.pdf"
    String quality_by_cycle_metrics = "~{output_bam_prefix}.quality_by_cycle_metrics"
    String quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    String quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
  }
}

# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
  input {
    String input_bam
    String input_bam_index
    String output_bam_prefix
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Boolean collect_gc_bias_metrics = true
    Int preemptible_tries
  }
  String flag_file = output_bam_prefix + ".done"  # finished flag
  Int memory_size = 7 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics

    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP
      
    touch "~{flag_file}"  # flag successful finish
  fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    preemptible: preemptible_tries
  }
  output {
    String alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    String gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    String gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    String gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
  }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
  input {
    String input_bam
    String input_bam_index
    String output_bam_prefix
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Boolean collect_gc_bias_metrics = true
    Int preemptible_tries
  }
  String flag_file = output_bam_prefix + ".done"  # finished flag
  Int memory_size = 7 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics \
      ~{output_bam_prefix}.insert_size_metrics \
      ~{output_bam_prefix}.insert_size_histogram.pdf

    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m  -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=QualityScoreDistribution \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY
      
    touch "~{flag_file}"  # flag successful finish
  fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    preemptible: preemptible_tries
  }
  output {
    String alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    String bait_bias_detail_metrics = "~{output_bam_prefix}.bait_bias_detail_metrics"
    String bait_bias_summary_metrics = "~{output_bam_prefix}.bait_bias_summary_metrics"
    String gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    String gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    String gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
    String insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    String insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    String pre_adapter_detail_metrics = "~{output_bam_prefix}.pre_adapter_detail_metrics"
    String pre_adapter_summary_metrics = "~{output_bam_prefix}.pre_adapter_summary_metrics"
    String quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    String quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
    String error_summary_metrics = "~{output_bam_prefix}.error_summary_metrics"
  }
}

task ConvertSequencingArtifactToOxoG {
  input {
    String pre_adapter_detail_metrics
    String bait_bias_detail_metrics
    String base_name
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Int preemptible_tries
    Int memory_multiplier = 1
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  # Int disk_size = ceil(size(pre_adapter_detail_metrics, "GiB") + size(bait_bias_detail_metrics, "GiB") + ref_size) + 20

  Int memory_size = ceil(5 * memory_multiplier)
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
    input_base=$(dirname ~{pre_adapter_detail_metrics})/~{base_name}
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -Dpicard.useLegacyParser=false \
      -jar /usr/picard/picard.jar \
      ConvertSequencingArtifactToOxoG \
      --INPUT_BASE $input_base \
      --OUTPUT_BASE ~{base_name} \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    preemptible: preemptible_tries
  }
  output {
    File oxog_metrics = "~{base_name}.oxog_metrics"
  }
}

# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  input {
    Array[String] input_bams
    Array[String] input_bam_indexes
    String haplotype_database_file
    String metrics_filename
    # Float total_input_size
    Int preemptible_tries
    Float lod_threshold
    String cross_check_by
  }

  String flag_file = metrics_filename + ".done"  # finished flag
  Int memory_size = 5
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command <<<
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{min_java_mem}m -Xmx~{java_memory_size}m \
      -jar /usr/picard/picard.jar \
      CrosscheckFingerprints \
      OUTPUT=~{metrics_filename} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      EXPECT_ALL_GROUPS_TO_MATCH=true \
      INPUT=~{sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=~{lod_threshold} \
      CROSSCHECK_BY=~{cross_check_by}

    touch "~{flag_file}" # flag successful finish
  fi
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }
  output {
    String cross_check_fingerprints_metrics = metrics_filename
  }
}

# Check that the fingerprint of the sample BAM matches the sample array
task CheckFingerprint {
  input {
    String input_bam
    String input_bam_index
    String output_basename
    String haplotype_database_file
    String? genotypes
    String? genotypes_index
    String sample
    Int preemptible_tries
  }
  # Picard has different behavior depending on whether or not the OUTPUT parameter ends with a '.', so we are explicitly
  #   passing in where we want the two metrics files to go to avoid any potential confusion.
  String summary_metrics_location = "~{output_basename}.fingerprinting_summary_metrics"
  String detail_metrics_location = "~{output_basename}.fingerprinting_detail_metrics"
  
  String flag_file = summary_metrics_location + ".done"
  
  Int memory_size = 5
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command <<<
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{min_java_mem}m -Xmx~{java_memory_size}m  \
      -jar /usr/picard/picard.jar \
      CheckFingerprint \
      INPUT=~{input_bam} \
      SUMMARY_OUTPUT=~{summary_metrics_location} \
      DETAIL_OUTPUT=~{detail_metrics_location} \
      GENOTYPES=~{genotypes} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      SAMPLE_ALIAS="~{sample}" \
      IGNORE_READ_GROUPS=true

    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  >>>
 runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
  }
  output {
    String summary_metrics = summary_metrics_location
    String detail_metrics = detail_metrics_location
  }
}

task CheckPreValidation {
  input {
    String duplication_metrics
    String chimerism_metrics
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    Int preemptible_tries
  }

  command <<<
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

python <<CODE

import csv
with open('duplication.csv') as dupfile:
  reader = csv.DictReader(dupfile, delimiter='\t')
  for row in reader:
    with open("duplication_value.txt","w") as file:
      file.write(row['PERCENT_DUPLICATION'])
      file.close()

with open('chimerism.csv') as chimfile:
  reader = csv.DictReader(chimfile, delimiter='\t')
  for row in reader:
    with open("chimerism_value.txt","w") as file:
      file.write(row['PCT_CHIMERAS'])
      file.close()

CODE

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "2 GiB"
  }
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}

task ValidateSamFile {
  input {
    String input_bam
    String? input_bam_index
    String report_filename
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Int? max_output
    Array[String]? ignore
    Boolean? is_outlier_data
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }

  String flag_file = report_filename + ".done"

  Int memory_size = ceil(7 * memory_multiplier)
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m  -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
      
    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }
  output {
    String report = "~{report_filename}"
  }
}

# Note these tasks will break if the read lengths in the bam are greater than 250.
task CollectWgsMetrics {
  input {
    String input_bam
    String input_bam_index
    String metrics_filename
    String wgs_coverage_interval_list
    String ref_fasta
    String ref_fasta_index
    Int read_length
    Int preemptible_tries
  }
  String flag_file = metrics_filename + ".done"
  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
      
    touch "~{flag_file}" # flag successful finish
  fi # test file exists

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
  }
  output {
    String metrics = "~{metrics_filename}"
  }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
  input {
    String input_bam
    String input_bam_index
    String metrics_filename
    String wgs_coverage_interval_list
    String ref_fasta
    String ref_fasta_index
    Int read_length
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }

  String flag_file = metrics_filename + ".done"
  Int memory_size = ceil(7 * memory_multiplier)  
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectRawWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
      
    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    cpu: "1"
  }
  output {
    String metrics = "~{metrics_filename}"
  }
}

task CollectHsMetrics {
  input {
    String input_bam
    String input_bam_index
    String ref_fasta
    String ref_fasta_index
    String metrics_filename
    String target_interval_list
    String bait_interval_list
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }


  Int memory_size = 10
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  String flag_file = metrics_filename + ".done"

  # There are probably more metrics we want to generate with this tool
  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectHsMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      VALIDATION_STRINGENCY=SILENT \
      TARGET_INTERVALS=~{target_interval_list} \
      BAIT_INTERVALS=~{bait_interval_list} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY \
      OUTPUT=~{metrics_filename}
      
    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }

  output {
    String metrics = metrics_filename
  }
}

# Generate a checksum per readgroup
task CalculateReadGroupChecksum {
  input {
    String input_bam
    String input_bam_index
    String read_group_md5_filename
    Int preemptible_tries
  }

  String flag_file = read_group_md5_filename + ".done"
  Int memory_size = 3
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m  -jar /usr/picard/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=~{input_bam} \
      OUTPUT=~{read_group_md5_filename}

    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }
  output {
    String md5_file = "~{read_group_md5_filename}"
  }
}

# Validate a (g)VCF with -gvcf specific validation
task ValidateVCF {
  input {
    String input_vcf
    String input_vcf_index
    String ref_fasta
    String ref_fasta_index
    String ref_dict
    String dbsnp_vcf
    String dbsnp_vcf_index
    String calling_interval_list
    Int preemptible_tries
    Boolean is_gvcf = true
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  String flag_file = input_vcf + ".validate.done"
  Int memory_size = 7 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    gatk --java-options " -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m "  \
      ValidateVariants \
      -V ~{input_vcf} \
      -R ~{ref_fasta} \
      -L ~{calling_interval_list} \
      ~{true="-gvcf" false="" is_gvcf} \
      --validation-type-to-exclude ALLELES \
      --dbsnp ~{dbsnp_vcf}
    touch ~{flag_file}
  fi # test file exists
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }
}

# Collect variant calling metrics from GVCF output
task CollectVariantCallingMetrics {
  input {
    String input_vcf
    String input_vcf_index
    String metrics_basename
    String dbsnp_vcf
    String dbsnp_vcf_index
    String ref_dict
    String evaluation_interval_list
    Boolean is_gvcf = true
    Int preemptible_tries
  }

  String flag_file = metrics_basename + ".variant_calling_qc.done"
  Int memory_size = 5 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=~{input_vcf} \
      OUTPUT=~{metrics_basename} \
      DBSNP=~{dbsnp_vcf} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS=~{evaluation_interval_list} \
      ~{true="GVCF_INPUT=true" false="" is_gvcf}
      
    touch ~{flag_file}
  fi # test file exists
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
  }
  output {
    String summary_metrics = "~{metrics_basename}.variant_calling_summary_metrics"
    String detail_metrics = "~{metrics_basename}.variant_calling_detail_metrics"
  }
}
