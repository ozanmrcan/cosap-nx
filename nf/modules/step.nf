/*
 * Placeholder process used for every stage while the graph wiring is tested.
 * Writes a small text file instead of running a real tool.
 */

def EXTENSIONS = [
    trim: 'fastq', mapping: 'bam', mdup: 'bam', calibrate: 'bam',
    variant_calling: 'vcf', msi: 'txt', genefusion: 'json'
]

process STEP {
    tag "${meta.stage}/${meta.id}"
    cpus   { meta.cpus ?: params.defaults.cpus }
    memory { meta.memory ?: params.defaults.memory }
    publishDir "${params.outdir}/${meta.stage}", mode: 'copy'

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("${meta.id}.${EXTENSIONS[meta.stage]}")

    script:
    """
    echo "${meta.library ?: meta.stage} on: ${inputs}" > ${meta.id}.${EXTENSIONS[meta.stage]}
    """
}
