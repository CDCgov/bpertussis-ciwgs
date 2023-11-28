process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    tuple val(meta), path("*.R{1,2}.fastq.gz")     , optional:true, emit: fastq
    tuple val(meta), path("*_interleaved.fastq.gz"), optional:true, emit: interleaved
    tuple val(meta), path("*_singleton.fastq.gz")  , optional:true, emit: singleton
    tuple val(meta), path("*_other.fastq.gz")      , optional:true, emit: other
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: ''
    def filt = task.ext.filt ?: ''
    def output = ( interleave && ! meta.single_end ) ? "> ${prefix}_interleaved.fastq.gz" :
        meta.single_end ? "-1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz" :
        "-1 ${prefix}.${filt}R1.fastq.gz -2 ${prefix}.${filt}R2.fastq.gz" // -s ${prefix}_singleton.fastq.gz"
    """
    samtools sort -n --threads ${task.cpus-1} $input -o ${prefix}.${suffix}namesort.bam
    samtools \\
        fastq \\
        $args \\
        --threads ${task.cpus-1} \\
        ${prefix}.${suffix}namesort.bam \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
//        -0 ${prefix}_other.fastq.gz \\
