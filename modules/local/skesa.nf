process SKESA {
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::skesa=2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa:2.4.0--he1c1bb9_0' :
        'quay.io/biocontainers/skesa:2.4.0--he1c1bb9_0' }"

    input:
    tuple val(meta) , path(reads)

    output:
    tuple val(meta), path("*.fasta")  , emit: fasta
    tuple val(meta), path("*.log")    , emit: log
    path "versions.yml"               , emit: versions

    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def reads_args = ""
    if (meta.single_end) {
        reads_args = "${reads}"
    } else {
        reads_args = "${reads[0]},${reads[1]}"
    }

    """
    skesa \\
        $args \\
        --cores ${task.cpus} \\
        --contigs_out ${prefix}.skesa.fasta \\
        --reads $reads_args \\
        2> ${prefix}.skesa.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skesa: \$(skesa --version 2>&1 | grep SKESA | sed 's/^.*SKESA //; s/ .*\$//')
    END_VERSIONS
    """
}


