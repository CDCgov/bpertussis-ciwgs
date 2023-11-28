/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBpertussisciwgs.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.input) { ch_samplesheet = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
ch_hs_index = params.hsref

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { PREPARE_BED } from '../subworkflows/local/prepare_bed'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                            } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { BOWTIE2_ALIGN as HS_ALIGN         } from '../modules/nf-core/modules/nf-core/bowtie2/align/main_hs'
include { BOWTIE2_ALIGN as BP_ALIGN         } from '../modules/nf-core/modules/nf-core/bowtie2/align/main_bp'
include { SAMTOOLS_FLAGSTAT as HS_FLAGSTAT  } from '../modules/nf-core/modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as BP_FLAGSTAT  } from '../modules/nf-core/modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_SORT                     } from '../modules/nf-core/modules/nf-core/samtools/sort/main'
include { SAMTOOLS_DEPTH                    } from '../modules/nf-core/modules/nf-core/samtools/depth/main'
include { SAMTOOLS_FASTQ as HS_FASTQ        } from '../modules/nf-core/modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as BP_FASTQ        } from '../modules/nf-core/modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_BEDCOV                   } from '../modules/local/samtools_bedcov'
include { SKESA                             } from '../modules/local/skesa'
include { QUAST                             } from '../modules/nf-core/modules/nf-core/quast/main'
include { SUMMARIZE_QC                      } from '../modules/local/summarize'
include { DUMPSOFTWAREVERSIONS              } from '../modules/local/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BPERTUSSISCIWGS {

    ch_software_versions = Channel.empty()

    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )

    // MODULE: Run FastQC
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))

    // SUBWORKFLOW: Prepare reference
    PREPARE_GENOME ()
    ch_bt_index = PREPARE_GENOME.out.bt_index
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.versions)

    // SUBWORKFLOW: Prepare MLST bed file if allele refs provided.
    if ( params.mlst ){
        PREPARE_BED ()
        ch_mlst_bed = PREPARE_BED.out.bedout
        ch_software_versions = ch_software_versions.mix(PREPARE_BED.out.versions)
    }

// The main CIWGS read filtering workflow starts here
// Step-1: Subtractive read mapping to human genome
    // MODULE: Run bowtie2-align to human reference
    HS_ALIGN (
        INPUT_CHECK.out.reads,
        ch_hs_index,
        false,
        true
    )
    ch_hs_unaln = HS_ALIGN.out.bam
    ch_hs_flagstat = HS_ALIGN.out.bam
    ch_software_versions = ch_software_versions.mix(HS_ALIGN.out.versions)

    // MODULE: Flagstat
    HS_FLAGSTAT (
        ch_hs_flagstat
    )
    ch_software_versions = ch_software_versions.mix(HS_FLAGSTAT.out.versions)

    // MODULE: Run samtools fastq, export non-human reads
    HS_FASTQ (
        ch_hs_unaln,
        false
    )
    ch_nonhs_reads = HS_FASTQ.out.fastq
    ch_software_versions = ch_software_versions.mix(HS_FASTQ.out.versions)

// Step-2: Positive filter read mapping to B. pertussis genome
    // MODULE: Run bowtie2 align to bp reference // 20230411 - could probably streamline this custom module copy
    BP_ALIGN (
        ch_nonhs_reads,
        ch_bt_index,
        false,
        true
    )
    ch_bp_aln = BP_ALIGN.out.bam
    ch_bp_flagstat = BP_ALIGN.out.bam
    ch_software_versions = ch_software_versions.mix(BP_ALIGN.out.versions)

    // MODULE: Flagstat
    BP_FLAGSTAT (
        ch_bp_flagstat
    )
    ch_software_versions = ch_software_versions.mix(BP_FLAGSTAT.out.versions)

    // MODULE: Run samtools sort (and index)
    SAMTOOLS_SORT (
        ch_bp_aln
    )
    ch_bp_sorted = SAMTOOLS_SORT.out.bam
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_SORT.out.versions)

    // MODULE: Run samtools fastq, export bp-filtered reads (first namesort bam file)
    BP_FASTQ (
        ch_bp_sorted,
        false
    )
    ch_bp_reads = BP_FASTQ.out.fastq
    ch_software_versions = ch_software_versions.mix(BP_FASTQ.out.versions)
    //can i add simple process to rename fastq files?

// Step-3: Evaluate target genome recovery
    // MODULE: Run samtools depth
    SAMTOOLS_DEPTH (
        ch_bp_sorted
    )
    ch_bp_depth = SAMTOOLS_DEPTH.out.tsv
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_DEPTH.out.versions)

    // If wgMLST/cgMLST/MLST allele fasta provided
    if ( params.mlst ) {
        
        SAMTOOLS_BEDCOV (
            ch_bp_sorted,
            ch_mlst_bed
        )
        ch_software_versions = ch_software_versions.mix(SAMTOOLS_BEDCOV.out.versions)

    }

    // **OPTIONAL**
    if ( params.skesa == true ) {
    // MODULE: De novo assembly with skesa
        SKESA (
            ch_bp_reads
        )
        ch_bp_contigs = SKESA.out.fasta
        ch_software_versions = ch_software_versions.mix(SKESA.out.versions)

        // MODULE: QUAST assembly QC
        QUAST (
            SKESA.out.fasta.collect{it[1]}.ifEmpty([]),
            params.fasta
        )
    }

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBpertussisciwgs.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

    //
    // MODULE: capture and summarize QC 
    //
    ch_summarizeqc_files = Channel.empty()
    ch_summarizeqc_files = ch_summarizeqc_files.mix(MULTIQC.out.data.collect().ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(SAMTOOLS_DEPTH.out.tsv.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(HS_FLAGSTAT.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(BP_FLAGSTAT.out.flagstat.collect{it[1]}.ifEmpty([]))
    if ( params.mlst ) {
        ch_summarizeqc_files = ch_summarizeqc_files.mix(SAMTOOLS_BEDCOV.out.tsv.collect{it[1]}.ifEmpty([]))
    }
    if ( params.skesa == true ){
        ch_summarizeqc_files = ch_summarizeqc_files.mix(QUAST.out.tsv.collect().ifEmpty([]))
    }
    SUMMARIZE_QC (
        ch_summarizeqc_files.collect(),
        ch_samplesheet
    )

// Step-4: flotsam and jetsam
    // MODULE: Pipeline reporting
    //
    DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )


}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
