process DECIMATE {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/decimated", mode: 'copy', pattern: "${meta.pulsar}_${meta.utc}_zap*t.ar"


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:latest' }"

    input:
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
