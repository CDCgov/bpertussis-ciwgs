//
// Generate reference genome related files for analysis
//

include { BLASTN_PAIRWISE  } from '../../modules/nf-core/modules/nf-core/blast/blastn/main_pairwise'


workflow PREPARE_BED {
    main:
    
    //
    // generate bed file for mlst locus coords in bp ref genome
    // 

    ch_versions = Channel.empty()
    ch_fasta = file(params.fasta)
    ch_mlst = file(params.mlst)

    BLASTN_PAIRWISE ( 
        ch_mlst,
        ch_fasta
    )
    
    ch_mlst_bls = BLASTN_PAIRWISE.out.bls
    ch_mlst_bed = BLASTN_PAIRWISE.out.bed
    ch_versions = ch_versions.mix(BLASTN_PAIRWISE.out.versions)

    emit:

    blastout    = ch_mlst_bls
    bedout      = ch_mlst_bed
    versions    = ch_versions         // channel: [ versions.yml ]

}
