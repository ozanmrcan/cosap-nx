#!/usr/bin/env nextflow

/*
 * COSAP-NX
 *
 * Reads a pipeline JSON where every stage lists named steps, and every step
 * names the step(s) it reads from as "stage/name". Stages run in a fixed order,
 * each stage can hold any number of steps.
 *
 * Tools with a real process: bwa, deepvariant, haplotypecaller.
 * Everything else uses the placeholder STEP process and only runs with -stub-run.
 */

nextflow.enable.dsl = 2

include { BWA             } from './modules/bwa'
include { DEEPVARIANT     } from './modules/deepvariant'
include { HAPLOTYPECALLER } from './modules/haplotypecaller'

include { STEP as TRIM          } from './modules/step'
include { STEP as MAPPING_STUB  } from './modules/step'
include { STEP as MDUP          } from './modules/step'
include { STEP as CALIBRATE     } from './modules/step'
include { STEP as CALLING_STUB  } from './modules/step'
include { STEP as MSI           } from './modules/step'
include { STEP as GENEFUSION    } from './modules/step'

// Inputs of a step by role: [input: 'trim/normal'] or [normal: '...', tumor: '...']
def stepInputs(spec) {
    spec.subMap(['input', 'normal', 'tumor']).findAll { it.value }
}

def stageSteps(stage) {
    params.containsKey(stage) ? params[stage] : [:]
}

/*
 * Finds the steps of a stage that read from the outputs seen so far.
 * Outputs come in as [ref, meta, files]. A step with two inputs (normal + tumor)
 * waits until both have arrived. Emits [meta, files] ready for the process,
 * files ordered normal first, then tumor.
 */
def link(stage, done) {
    def steps = stageSteps(stage)

    done
        .flatMap { ref, meta, files ->
            def hits = []
            steps.each { name, spec ->
                def ins = stepInputs(spec)
                ins.each { role, r ->
                    if (r == ref) {
                        hits << [groupKey(name, ins.size()), role, meta, files]
                    }
                }
            }
            hits
        }
        .groupTuple()
        .map { key, roles, metas, files ->
            def name = key.toString()
            def spec = steps[name]
            def up = [roles, metas].transpose().collectEntries()
            def upFiles = [roles, files].transpose().collectEntries()

            def meta = [
                id     : name,
                stage  : stage,
                library: spec.library,
                params : spec.params ?: [:],
                cpus   : spec.cpus,
                memory : spec.memory,
                spark  : spec.spark,
                gvcf   : spec.gvcf ?: false,
                // labels carried forward, the comparator needs mapper + caller
                sample : up.input?.sample ?: up.normal?.sample,
                normal : up.normal?.sample,
                tumor  : up.tumor?.sample,
                mapper : stage == 'mapping' ? spec.library : (up.input ?: up.normal ?: up.tumor).mapper,
                caller : stage == 'variant_calling' ? spec.library : null,
            ]
            def ordered = ['input', 'normal', 'tumor'].findAll { upFiles[it] }.collect { upFiles[it] }
            [meta, ordered.flatten()]
        }
}

// Process output [meta, file, file, ...] -> [ref, meta, files] so later stages can find it
def asRef(ch) {
    ch.map { row -> ["${row[0].stage}/${row[0].id}".toString(), row[0], row[1..-1].flatten()] }
}

// Input files of one entry in "inputs": FASTQ pair, or BAM + index
def inputFiles(name, s) {
    if (s.type == 'fastq') {
        return [file(s.r1, checkIfExists: !workflow.stubRun), file(s.r2, checkIfExists: !workflow.stubRun)]
    }
    if (s.type == 'bam') {
        def bai = [file("${s.bam}.bai"), file(s.bam.toString().replaceAll(/\.bam$/, '.bai'))].find { it.exists() }
        if (!bai && !workflow.stubRun) error "inputs/${name}: no index found for ${s.bam}"
        return [file(s.bam, checkIfExists: !workflow.stubRun)] + (bai ? [bai] : [])
    }
    error "inputs/${name}: unknown input type '${s.type}'"
}

workflow {
    def STAGES = ['trim', 'mapping', 'mdup', 'calibrate', 'variant_calling', 'msi', 'genefusion']
    // Tools with a real process, per stage
    def TOOLS = [mapping: ['bwa'], variant_calling: ['deepvariant', 'haplotypecaller']]

    // Check the graph before running anything
    def known = params.inputs.keySet().collect { "inputs/${it}".toString() }
    STAGES.each { stage -> stageSteps(stage).each { name, spec -> known << "${stage}/${name}".toString() } }

    STAGES.each { stage ->
        stageSteps(stage).each { name, spec ->
            stepInputs(spec).each { role, ref ->
                if (!(ref in known)) error "${stage}/${name}: ${role} '${ref}' is not a step in this pipeline"
            }
            if (!workflow.stubRun && !(spec.library in (TOOLS[stage] ?: []))) {
                error "${stage}/${name}: '${spec.library ?: stage}' is not supported yet (only runs with -stub-run)"
            }
        }
    }

    // Reference fasta plus everything next to it (.fai, .dict, bwa index)
    def fasta_path = file(params.references.fasta, checkIfExists: !workflow.stubRun)
    def fasta = Channel.value(fasta_path)
    def ref_index = Channel.value(
        files("${fasta_path}.*") + files("${fasta_path.parent}/${fasta_path.baseName}.dict")
    )

    def done = Channel.fromList(
        params.inputs.collect { name, s ->
            ["inputs/${name}".toString(), [sample: name], inputFiles(name, s)]
        }
    )

    done = done.mix(asRef(TRIM(link('trim', done))))

    def mapping = link('mapping', done).branch { meta, files ->
        bwa  : meta.library == 'bwa'
        other: true
    }
    def mapped = BWA(mapping.bwa, fasta, ref_index)
        .mix(MAPPING_STUB(mapping.other))
    done = done.mix(asRef(mapped))

    done = done.mix(asRef(MDUP(link('mdup', done))))
    done = done.mix(asRef(CALIBRATE(link('calibrate', done))))

    def calling = link('variant_calling', done).branch { meta, files ->
        deepvariant    : meta.library == 'deepvariant'
        haplotypecaller: meta.library == 'haplotypecaller'
        other          : true
    }
    DEEPVARIANT(calling.deepvariant, fasta, ref_index)
    def called = DEEPVARIANT.out[0]
        .mix(HAPLOTYPECALLER(calling.haplotypecaller, fasta, ref_index))
        .mix(CALLING_STUB(calling.other))
    done = done.mix(asRef(called))

    MSI(link('msi', done))
    GENEFUSION(link('genefusion', done))
}
