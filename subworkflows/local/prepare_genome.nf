//
// Generate reference genome related files for analysis
//

include { BOWTIE2_BUILD } from '../../modules/nf-core/modules/nf-core/bowtie2/build/main'


workflow PREPARE_GENOME {
    main:
    
    //
    // generate index for bwa
    // 

    ch_versions = Channel.empty()
    ch_fasta = file(params.fasta)

    BOWTIE2_BUILD ( ch_fasta )
    
    ch_bt_index = BOWTIE2_BUILD.out.index
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    emit:
    
    fasta                = ch_fasta            // path: genome.fasta
    bt_index            = ch_bt_index        // path: bowtie2/index/
    versions             = ch_versions         // channel: [ versions.yml ]
}
