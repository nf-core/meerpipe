// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process DECIMATE {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/decimated", mode: 'copy', pattern: "${meta.pulsar}_${meta.utc}_zap*t.ar"


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(cleaned_archive)

    output:
    tuple val(meta), path("${meta.pulsar}_${meta.utc}_zap*t.ar")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    if [ "${params.chop_edge}" == "true" ]; then
        chop_edge_channels --band ${meta.band} ${cleaned_archive}
        clean_ar=${cleaned_archive.getName().replace("_zap", "_zap_chopped")}
    else
        clean_ar=${cleaned_archive}
    fi
    for nchan in ${meta.nchans.join(' ')}; do
        if ${params.use_max_nsub}; then
            nsub=\$( vap -c nsub  \$clean_ar | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
            # Calculate nsub to get desired TOA S/N
            nsubs=\$(calc_max_nsub --sn ${meta.snr} --nchan \${nchan} --duration ${meta.dur} --input_nsub \${nsub} --sn_desired ${params.tos_sn})
        else
            nsubs="1"
        fi

        # Make a max_nsub decimation and a time scrunched decimation
        for nsub in \$nsubs; do
            # Make full stokes and/or polarisation scrunched
            for stokes in ${meta.npols.join(' ')}; do
                if [ \${stokes} -eq 1 ]; then
                    # Polarisation scrunch option
                    stokes_op="-p"
                else
                    stokes_op=""
                fi

                echo "Decimate nsub=\${nsub}  nchan=\${nchan} stokes=\${stokes}"
                pam --setnsub \${nsub} --setnchn \${nchan} -S \${stokes_op} -e \${nchan}ch\${stokes}p\${nsub}t.ar \${clean_ar}
            done
        done
    done
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${meta.pulsar}_${meta.utc}_zap.1ch1p1t.ar
    """
}
