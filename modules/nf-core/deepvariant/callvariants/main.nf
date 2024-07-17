
process DEEPVARIANT_CALLVARIANTS {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.5.0"

    input:
    tuple val(make_examples_id), val(meta), path(make_examples_tfrecords)
    val model_type

    output:
    tuple val(make_examples_id), val(meta), path("${prefix}.call.tfrecord.gz"),     emit: call_variants_tfrecords
    path "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def matcher = make_examples_tfrecords[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    if (!matcher.matches()) {
        throw new IllegalArgumentException("tfrecord baseName '" + make_examples_tfrecords[0].baseName + "' doesn't match the expected pattern")
    }
    def examples_tfrecord_name = matcher[0][1]
    def shardCount = matcher[0][2]
    // Reconstruct the logical name - ${tfrecord_name}.examples.tfrecord@${task.cpus}.gz
    def examples_tfrecords_logical_name = "${examples_tfrecord_name}@${shardCount}.gz"
    
    def model_type_clean = model_type.toString().replaceAll("[^A-Za-z0-9_-]", "")
    """
    /opt/deepvariant/bin/call_variants \\
        --outfile ${prefix}.call.tfrecord.gz \\
        --examples ${examples_tfrecords_logical_name} \\
        --checkpoint "/opt/models/${model_type_clean}/model.ckpt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_callvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.call.tfrecord.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_callvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

}
