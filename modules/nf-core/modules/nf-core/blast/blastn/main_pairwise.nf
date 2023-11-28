process BLASTN_PAIRWISE {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::blast=2.13.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0' :
        'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    //tuple val(meta), path(fasta)
    //path  db
    path(query)
    path(fasta)

    output:
    path('*.mlst.bls'),   emit: bls
    path('*.mlst.bed'),   emit: bed
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    blastn \\
        -subject $fasta \\
        -query $query \\
        $args \\
        -outfmt '6 sseqid sstart send qseqid' \\
        -out ${prefix}.mlst.bls
    cat ${prefix}.mlst.bls | awk '{OFS="\t"; if(\$3<\$2) {print \$1,\$3,\$2,\$4} else {print \$0} }' | sort -n -k2 > ${prefix}.mlst.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
