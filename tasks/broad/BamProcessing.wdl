version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for BAM file processing of human whole-genome or exome sequencing data.
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

# Sort BAM file by coordinate order
task SortSam {
  input {
    String input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
  }
  
  String flag_file = output_bam_basename + ".done"  # finished flag

  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Int memory_size = 5 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
    set -e
    if [ -f "~{flag_file}" ]; then
      echo "SKIP - file ~{flag_file} already exists."
    else
    java -Dsamjdk.compression_level=~{compression_level} -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

      touch "~{flag_file}" # flag successful finish
    fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    cpu: "1"
    memory: "~{memory_size} GiB"
    preemptible: preemptible_tries
  }
  output {
    String output_bam = "~{output_bam_basename}.bam"
    String output_bam_index = "~{output_bam_basename}.bai"
    String output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Sort BAM file by coordinate order -- using Spark!
task SortSamSpark {
  input {
    String input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier

  command {
    set -e

    gatk --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms100g -Xmx100g" \
      SortSamSpark \
      -I ~{input_bam} \
      -O ~{output_bam_basename}.bam \
      -- --conf spark.local.dir=. --spark-master 'local[16]' --conf 'spark.kryo.referenceTracking=false'

    samtools index ~{output_bam_basename}.bam ~{output_bam_basename}.bai
  }
  runtime {
    docker: gatk_docker
    bootDiskSizeGb: "15"
    cpu: "16"
    memory: "102 GiB"
  }
  output {
    String output_bam = "~{output_bam_basename}.bam"
    String output_bam_index = "~{output_bam_basename}.bai"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[String] input_bams
    String output_bam_basename
    String metrics_filename
    # Float total_input_size
    Int compression_level
    # Int preemptible_tries
    String tmp_dir

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk = 20
  }

  String flag_file = output_bam_basename + ".done"  # finished flag

  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom

  Float memory_size = 10 * memory_multiplier
  Int java_memory_size = ceil(memory_size - 1)
  Int min_mem_size = ceil(java_memory_size / 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    set -e
    if [ -f "~{flag_file}" ]; then
      echo "SKIP - file ~{flag_file} already exists."
    else
    java -Dsamjdk.compression_level=~{compression_level} -server -XX:+UseG1GC -Xms~{min_mem_size}g -Xmx~{java_memory_size}g -server \
      -jar /usr/picard/picard.jar \
      MarkDuplicates \
      TMP_DIR=~{tmp_dir} \
      MAX_RECORDS_IN_RAM=4000000 \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false

      touch "~{flag_file}" # flag successful finish
    fi
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    cpu: "1"
  }
  output {
    String output_bam = "~{output_bam_basename}.bam"
    String duplicate_metrics = "~{metrics_filename}"
  }
}

task MarkDuplicatesSpark {
  input {
    Array[String] input_bams
    String output_bam_basename
    String metrics_filename
    Float total_input_size
    Int compression_level
    Int preemptible_tries

    String? read_name_regex
    Int memory_multiplier = 3
    Int cpu_size = 6
  }

  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving 2.5 as wiggleroom
  Float md_disk_multiplier = 2.5
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + 20

  Int memory_size = ceil(16 * memory_multiplier)
  Int java_memory_size = (memory_size - 6)

  String output_bam_location = "~{output_bam_basename}.bam"

  # Removed options ASSUME_SORT_ORDER, CLEAR_DT, and ADD_PG_TAG_TO_READS as it seems like they are a) not implemented
  #   in MarkDuplicatesSpark, and/or b) are set to "false" aka "don't do" anyhow.
  # MarkDuplicatesSpark requires PAPIv2
  command <<<
    set -e
    export GATK_LOCAL_JAR=/root/gatk.jar
    gatk --java-options "-Dsamjdk.compression_level=~{compression_level} -Xmx~{java_memory_size}g" \
      MarkDuplicatesSpark \
      --input ~{sep=' --input ' input_bams} \
      --output ~{output_bam_location} \
      --metrics-file ~{metrics_filename} \
      --read-validation-stringency SILENT \
      ~{"--read-name-regex " + read_name_regex} \
      --optical-duplicate-pixel-distance 2500 \
      --treat-unsorted-as-querygroup-ordered \
      --create-output-bam-index false \
      -- --conf spark.local.dir=/mnt/tmp --spark-master 'local[16]' --conf 'spark.kryo.referenceTracking=false'
  >>>

  runtime {
    docker: "jamesemery/gatknightly:gatkMasterSnapshot44ca2e9e84a"
    cpu: cpu_size
    memory: "~{memory_size} GiB"
  }

  output {
    String output_bam = output_bam_location
    String duplicate_metrics = metrics_filename
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    String input_bam
    String recalibration_report_filename
    Array[String] sequence_group_interval
    String dbsnp_vcf
    String dbsnp_vcf_index
    Array[String] known_indels_sites_vcfs
    Array[String] known_indels_sites_indices
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  String recalibration_report_seq_group_filename = recalibration_report_filename + "." + sequence_group_interval[0]
  String flag_file = recalibration_report_seq_group_filename + ".done"  # finished flag

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    set -e
    echo "~{sep=';' sequence_group_interval}"
    if [ -f "~{flag_file}" ]; then
      echo "SKIP - file ~{flag_file} already exists."
    else
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g -Xmx7g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O "~{recalibration_report_seq_group_filename}" \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}

      touch "~{flag_file}" # flag successful finish
    fi
  }
  runtime {
    docker: gatk_docker
    memory: "8 GiB"
    cpu: "1"
  }
  output {
    String recalibration_report = "~{recalibration_report_seq_group_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    String input_bam
    String output_bam_basename
    String recalibration_report
    Array[String] sequence_group_interval
    String ref_dict
    String ref_fasta
    String ref_fasta_index
    Int compression_level
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
    Int memory_multiplier = 1
    Int additional_disk = 20
    Boolean bin_base_qualities = true
    Boolean somatic = false
  }

  String output_bam_basename_seq_group = output_bam_basename + "." + sequence_group_interval[0]
  String flag_file = output_bam_basename_seq_group + ".done"  # finished flag


  Int java_memory_size = ceil(3500 * memory_multiplier)
  Int memory_size = java_memory_size + 1

  Boolean bin_somatic_base_qualities = bin_base_qualities && somatic

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx~{memory_size}m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename_seq_group}.bam \
      -bqsr ~{recalibration_report} \
      ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 40' false='' bin_somatic_base_qualities} \
      ~{true='--static-quantized-quals 50' false='' bin_somatic_base_qualities} \
      -L ~{sep=" -L " sequence_group_interval}

    touch "~{flag_file}" # flag successful finish
  fi # test file exists
  }
  runtime {
    docker: gatk_docker
    memory: "~{memory_size} MiB"
    cpu: "1"
  }
  output {
    String recalibrated_bam = "~{output_bam_basename_seq_group}.bam"
    String recalibrated_bam_checksum = "~{output_bam_basename_seq_group}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  input {
    Array[String] input_bqsr_reports
    String output_report_filename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  String flag_file = output_report_filename + ".done"  # finished flag

  command {
    set -e
    if [ -f "~{flag_file}" ]; then
      echo "SKIP - file ~{flag_file} already exists."
    else
    gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}

      touch "~{flag_file}" # flag successful finish
    fi # test file exists
    }
  runtime {
    docker: gatk_docker
    memory: "3500 MiB"
    cpu: "1"
  }
  output {
    String output_bqsr_report = "~{output_report_filename}"
  }
}

# Combine multiple *sorted* BAM files
task GatherSortedBamFiles {
  input {
    Array[String] input_bams
    String output_bam_basename
    # Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  String flag_file = output_bam_basename + ".done"  # finished flag
  
  Int memory_size = 3 
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    java -Dsamjdk.compression_level=~{compression_level} -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m  -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true

    touch "~{flag_file}" # flag successful finish
  fi # test file exists
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    cpu: "1"
  }
  output {
    String output_bam = "~{output_bam_basename}.bam"
    String output_bam_index = "~{output_bam_basename}.bai"
    String output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple *unsorted* BAM files
# Note that if/when WDL supports optional outputs, we should merge this task with the sorted version
task GatherUnsortedBamFiles {
  input {
    Array[String] input_bams
    String output_bam_basename
    # Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  Int memory_size = 3
  Int java_memory_size = (memory_size - 1) * 1000
  Int min_java_mem = ceil(java_memory_size / 2)

  command {
    set -e
    java -Dsamjdk.compression_level=~{compression_level} -XX:+UseG1GC -Xms~{min_java_mem}m -Xmx~{java_memory_size}m -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "~{memory_size} GiB"
    cpu: "1"
  }
  output {
    String output_bam = "~{output_bam_basename}.bam"
  }
}

task GenerateSubsettedContaminationResources {
  input {
    String bait_set_name
    String target_interval_list
    String contamination_sites_ud
    String contamination_sites_bed
    String contamination_sites_mu
    Int preemptible_tries
  }

  String output_ud = bait_set_name + "." + basename(contamination_sites_ud)
  String output_bed = bait_set_name + "." + basename(contamination_sites_bed)
  String output_mu = bait_set_name + "." + basename(contamination_sites_mu)
  String target_overlap_counts = "target_overlap_counts.txt"

  String flag_file = output_bed + ".done"  # finished flag

  command <<<
  set -e
  if [ -f "~{flag_file}" ]; then
    echo "SKIP - file ~{flag_file} already exists."
  else
    set -e -o pipefail

    grep -vE "^@" ~{target_interval_list} |
       awk -v OFS='\t' '$2=$2-1' |
       /app/bedtools intersect -c -a ~{contamination_sites_bed} -b - |
       cut -f6 > ~{target_overlap_counts}

    function restrict_to_overlaps() {
        # print lines from whole-genome file from loci with non-zero overlap
        # with target intervals
        WGS_FILE=$1
        EXOME_FILE=$2
        paste ~{target_overlap_counts} $WGS_FILE |
            grep -Ev "^0" |
            cut -f 2- > $EXOME_FILE
        echo "Generated $EXOME_FILE"
    }

    restrict_to_overlaps ~{contamination_sites_ud} ~{output_ud}
    restrict_to_overlaps ~{contamination_sites_bed} ~{output_bed}
    restrict_to_overlaps ~{contamination_sites_mu} ~{output_mu}

    touch "~{flag_file}" # flag successful finish
  fi
  >>>
  runtime {
    memory: "3.5 GiB"
    docker: "us.gcr.io/broad-gotc-prod/bedtools:2.27.1"
    cpu: "1"
  }
  output {
    String subsetted_contamination_ud = output_ud
    String subsetted_contamination_bed = output_bed
    String subsetted_contamination_mu = output_mu
  }
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  input {
    String input_bam
    String input_bam_index
    String contamination_sites_ud
    String contamination_sites_bed
    String contamination_sites_mu
    String ref_fasta
    String ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Float contamination_underestimation_factor
    Boolean disable_sanity_check = false
  }

  String flag_file = output_prefix + ".done"  # finished flag

  command <<<
    set -e
    if [ ! -f "~{flag_file}" ]; then
    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null
    fi
      # used to read from the selfSM file and calculate contamination, which gets printed out
python3 <<CODE
import csv
import sys
with open('~{output_prefix}.selfSM') as selfSM:
  reader = csv.DictReader(selfSM, delimiter='\t')
  i = 0
  for row in reader:
    if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
      # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
      # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
      # vcf and bam.
      sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
      sys.exit(1)
    print(float(row["FREEMIX"])/~{contamination_underestimation_factor})
    i = i + 1
    # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
    # and the results are not reliable.
    if i != 1:
      sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
      sys.exit(2)
CODE
    touch "~{flag_file}" # flag successful finish
  >>>
  runtime {
    memory: "7.5 GiB"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c1cba76e979904eb69c31520a0d7f5be63c72253-1553018888"
    cpu: "2"
  }
  output {
    String selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}
