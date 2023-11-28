process SUMMARIZE_QC {

    //conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.22.2.1' :
        'quay.io/biocontainers/perl:5.22.2.1' }"

    input:
    path summarize_qc_files
    path samplesheet

    output:
    path '*.tsv'       , emit: tsv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in ./bin/
    """
    cp multiqc_data/multiqc_general_stats.txt ./
    
    summarize_qc.pl \\
        -dir . \\
        -samplesheet $samplesheet \\
        -out sample_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(perl --version | grep -o "v5")
    END_VERSIONS
    """
}
