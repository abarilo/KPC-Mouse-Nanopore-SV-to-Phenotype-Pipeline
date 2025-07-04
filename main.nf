#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
  //
  // 1) Input channels
  //
  Channel.fromPath(params.fastqs)
         .ifEmpty { error "No FASTQ files found: ${params.fastqs}" }
         .set { fastq_ch }

  Channel.value(file(params.ref))
         .ifEmpty { error "Reference not found: ${params.ref}" }
         .set { ref_ch }

  //
  // 2) QC → NanoPlot
  //
  nanoPlot(fastq_ch)

  //
  // 3) Alignment → BAMs
  //
  def bam_ch = align(fastq_ch, ref_ch)

  //
  // 4) Statistics → stats
  //
  collect_stats(bam_ch)

  //
  // 5a) Per-sample SNF (for population mode)
  //
  def snf_ch = sniffles2(bam_ch, ref_ch)

  //
  // 5b) Per-sample mosaic VCF
  //
  def mosaic_ch = snifflesMosaic(bam_ch, ref_ch)

  //
  // 6) Combined multi-sample VCF
  //
  def snf_list = snf_ch.map { sample, snf -> snf }.collect()
  def pop_ch   = snifflesPopulation(snf_list, ref_ch)

  //
  // 7) Plot both mosaic and population VCFs
  //
  snifflesPlotSample(mosaic_ch)
  snifflesPlotPop(pop_ch)


  //
 // 8) Annotate per-sample mosaic VCFs with Ensembl VEP
  //
  mosaic_ch
    .map { sample, vcf -> tuple(sample, vcf) }
    | vepMosaic

  //
  // 9) Annotate merged multi-sample VCF with Ensembl VEP
  //
  pop_ch
    .map { _, vcf -> vcf }
    | vepPopulation
}


////////////////////////////////////////////////////////////////////////////////
// Processes
////////////////////////////////////////////////////////////////////////////////

process nanoPlot {
  tag { fastq.baseName }
  publishDir "${params.outdir}/nanoplot", mode: 'copy'

  input:
    path fastq

  output:
    path "**"

  script:
  """
  base=\$(basename ${fastq} .fastq)
  NanoPlot \\
    --fastq ${fastq} \\
    --loglength \\
    --outdir . \\
    --prefix \$base
  """
}

process align {
  tag { fastq.baseName }
  cache true
  publishDir "${params.outdir}/bams", mode: 'copy'

  input:
    path fastq
    path ref

  output:
    tuple val(fastq.baseName),
          path("${fastq.baseName}.sorted.bam"),
          path("${fastq.baseName}.sorted.bam.bai")

  script:
  """
  REF=\$(basename ${ref} .gz)
  gunzip -c ${ref} > \$REF
  minimap2 -t ${task.cpus} -ax map-ont \$REF ${fastq} \\
    | samtools view -b -F 4 - \\
    | samtools sort -@ ${task.cpus} -o ${fastq.baseName}.sorted.bam
  samtools index ${fastq.baseName}.sorted.bam
  """
}

process collect_stats {
  tag { sample }
  publishDir "${params.outdir}/stats", mode: 'copy'

  input:
    tuple val(sample), path(bam), path(bai)

  output:
    path "${sample}.flagstat.txt"
    path "${sample}.stats.txt"

  script:
  """
  samtools flagstat ${bam}   > ${sample}.flagstat.txt
  samtools stats    ${bam}   > ${sample}.stats.txt
  """
}

process sniffles2 {
  tag { sample }
  publishDir "${params.outdir}/sv/snf", mode: 'copy'

  input:
    tuple val(sample), path(bam), path(bai)
    path ref

  output:
    tuple val(sample), path("${sample}.snf")

  script:
  """
  sniffles \\
    --input     ${bam} \\
    --snf       ${sample}.snf \\
    --reference ${ref} \\
    --threads   ${task.cpus}
  """
}

process snifflesMosaic {
  tag { sample }
  publishDir "${params.outdir}/sv/mosaic", mode: 'copy'

  input:
    tuple val(sample), path(bam), path(bai)
    path ref

  output:
    tuple val(sample), path("${sample}.mosaic.vcf")

  script:
  """
  REF_FA=\$(basename ${ref} .gz)
  gunzip -c ${ref} > \$REF_FA

  sniffles \\
    --input     ${bam} \\
    --vcf       ${sample}.mosaic.vcf \\
    --mosaic \\
    --reference \$REF_FA \\
    --threads   ${task.cpus}
  """
}

process snifflesPopulation {
  publishDir "${params.outdir}/sv/population", mode: 'copy'

  input:
    path snf_list    // collected list of all .snf paths
    path ref

  output:
    tuple val('population'), path("multisample.vcf")

  script:
  """
  REF_FA=\$(basename ${ref} .gz)
  gunzip -c ${ref} > \$REF_FA

  sniffles \\
    --input     ${snf_list.join(' ')} \\
    --vcf       multisample.vcf \\
    --reference \$REF_FA \\
    --threads   ${task.cpus}
  """
}

process snifflesPlotSample {
  tag { sample }
  publishDir "${params.outdir}/svplots/mosaic", mode: 'copy'

  input:
    tuple val(sample), path(vcf)

  output:
    path "**"

  script:
  """
  python3 -m sniffles2_plot \\
    -i ${vcf} \\
    -o .
  """
}

process snifflesPlotPop {
  tag { sample }
  publishDir "${params.outdir}/svplots/population", mode: 'copy'

  input:
    tuple val(sample), path(vcf)

  output:
    path "**"

  script:
  """
  python3 -m sniffles2_plot \\
    -i ${vcf} \\
    -o .
  """
}

process vepMosaic {
  tag { sample }
  publishDir "${params.outdir}/sv/vep/mosaic", mode: 'copy'

  input:
    tuple val(sample), path(mosaicVcf)

  output:
    tuple val(sample), path("${sample}.mosaic.vep.tsv")

  script:
  """
  vep \
    --species       ${params.vepSpecies} \
    --assembly      ${params.genomeBuild} \
    --cache --offline \
    --dir_cache     /root/.vep \
    --fasta         /root/.vep/mus_musculus/114_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa.gz \
    --tab \
    --fields        Location,Allele,Consequence,Gene,Symbol,Feature,BIOTYPE \
    --force_overwrite \
    --input_file    ${mosaicVcf} \
    --output_file   ${sample}.mosaic.vep.tsv
  """
}

process vepPopulation {
  publishDir "${params.outdir}/sv/vep/population", mode: 'copy'

  input:
    path populationVcf

  output:
    path("multisample.vep.tsv")

  script:
  """
  vep \
    --species       ${params.vepSpecies} \
    --assembly      ${params.genomeBuild} \
    --cache --offline \
    --dir_cache     /root/.vep \
    --fasta         /root/.vep/mus_musculus/114_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa.gz \
    --tab \
    --fields        Location,Allele,Consequence,Gene,Symbol,Feature,BIOTYPE \
    --force_overwrite \
    --input_file    ${populationVcf} \
    --output_file   multisample.vep.tsv
  """
}

