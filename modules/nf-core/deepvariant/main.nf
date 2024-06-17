process MAKE_EXAMPLES {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.5.0"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path(intervals),
            val("make_examples.tfrecord@${task.cpus}.gz"), path("make_examples.tfrecord-*-of-*.gz"),
            val("gvcf.tfrecord@${task.cpus}.gz"), path("gvcf.tfrecord-*-of-*.gz"),      emit: make_examples_tfrecords

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def regions = intervals ? "--regions=${intervals}" : ""

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \\
        --mode=calling \\
        --ref=${fasta} \\
        --reads=${input} \\
        --examples "./make_examples.tfrecord@${task.cpus}.gz" \\
        --channels "insert_size" \\
        --gvcf "./gvcf.tfrecord@${task.cpus}.gz" \\
        ${regions} \\
        ${args} \\
        --task {}
    """
}

process CALL_VARIANTS {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.5.0"

    input:
    tuple val(meta), path(intervals), val(make_example_tfrecord_filename), path(make_examples_tfrecords),
                                            val(gvcf_tfrecords_filename), path(gvcf_tfrecords)
    val model_type

    output:
    tuple val(meta), path(intervals), path("call_variants_output.tfrecord.gz"), val(gvcf_tfrecords_filename), path(gvcf_tfrecords),  emit: call_variants_tfrecords
    path "versions.yml"                                                         ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def model_type_clean = model_type.replaceAll("[^A-Za-z0-9]", "")

    """
    /opt/deepvariant/bin/call_variants \\
        --outfile call_variants_output.tfrecord.gz \\
        --examples ${make_example_tfrecord_filename} \\
        --checkpoint "/opt/models/${model_type_clean}/model.ckpt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

process POSTPROCESS_VARIANTS {
    tag "$meta.id"
    // TODO: This should be changed to a larger type when DeepVariant is updated to version 1.6
    label 'process_single'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.5.0"

    input:
    tuple val(meta), path(intervals), path(variant_calls_tfrecord), val(gvcf_tfrecords_filename), path(gvcf_tfrecords)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/deepvariant/bin/postprocess_variants \\
        --ref=${fasta} \\
        --infile=${variant_calls_tfrecord} \\
        --outfile=${prefix}.vcf.gz \\
        --nonvariant_site_tfrecord_path ${gvcf_tfrecords_filename} \\
        --gvcf_outfile=${prefix}.g.vcf.gz \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi
    """
}

workflow DEEPVARIANT {
    take:
    ch_input            // channel: [ val(meta), path(input), path(index), path(intervals)]
    ch_fasta            // channel: [ val(meta2), path(fasta) ]
    ch_fai              // channel: [ val(meta3), path(fail) ]
    ch_gzi              // channel: [ val(meta4), path(gzi) ]
    ch_model_type       // channel: val("wgs" | "wes")

    main:

    MAKE_EXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi)
    CALL_VARIANTS(MAKE_EXAMPLES.out.make_examples_tfrecords, ch_model_type)
    POSTPROCESS_VARIANTS(CALL_VARIANTS.out.call_variants_tfrecords, ch_fasta, ch_fai, ch_gzi)

    emit:
    vcf         = POSTPROCESS_VARIANTS.out.vcf
    vcf_tbi     = POSTPROCESS_VARIANTS.out.vcf_tbi
    gvcf        = POSTPROCESS_VARIANTS.out.gvcf
    gvcf_tbi    = POSTPROCESS_VARIANTS.out.gvcf_tbi
    versions    = CALL_VARIANTS.out.versions
}
