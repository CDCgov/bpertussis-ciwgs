#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/bpertussisciwgs
========================================================================================
    Github : https://github.com/cdcgov/bpertussis-ciwgs
    Website: https://github.com/cdcgov/bpertussis-ciwgs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { BPERTUSSISCIWGS } from './workflows/bpertussisciwgs'

//
// WORKFLOW: Run main bpertussisciwgs analysis pipeline
//
workflow NFCORE_BPERTUSSISCIWGS {
    BPERTUSSISCIWGS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// 
//
workflow {
    NFCORE_BPERTUSSISCIWGS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
