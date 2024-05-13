process DECIMATE {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/decimated", mode: 'copy', pattern: "${meta.pulsar}_${meta.utc}_{raw,zap}*t.ar"


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:3.0.6' }"

    input:
    tuple val(meta), path(template), path(archive)

    output:
    tuple val(meta), path(template), path("${meta.pulsar}_${meta.utc}_{raw,zap}*t.ar")

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
        chop_edge_channels --band ${meta.band} ${archive}
        archive=${archive.getName().replace("_zap", "_zap_chopped").replace("_raw", "_raw_chopped")}
    else
        archive=${archive}
    fi
    nsubs_list="1"
    if ${params.use_all_nsub}; then
        # Use all nsubs, do not time scrunch
        nsub=\$( vap -c nsub  \$archive | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        nsubs_list="\${nsubs_list} all_\${nsub}"
    fi
    if ${params.use_mode_nsub}; then
        # Use the most common observation length as the sub integration length
        mode="${ Math.floor( meta.dur.toFloat() / meta.mode_dur.toFloat() ).toInteger() }"
        if [ \${mode} -eq "0" ]; then
            mode=1
        fi
        nsubs_list="\${nsubs_list} mode_\${mode}"
    fi

    for nchan in ${meta.nchans.join(' ')}; do
        if [[ ${params.use_max_nsub} && "${meta.snr}" != "None" ]]; then
            nsub=\$( vap -c nsub  \$archive | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
            # Calculate nsub to get desired TOA S/N
            nsubs="\$nsubs_list \$(calc_max_nsub --sn ${meta.snr} --nchan \${nchan} --duration ${meta.dur} --input_nsub \${nsub} --sn_desired ${params.tos_sn})"
        else
            nsubs="\$nsubs_list"
        fi

        # Make a max_nsub decimation and a time scrunched decimation
        for nsub in \$nsubs; do
            if [[ \${nsub} == "all"* && \${nchan} -ne 1 ]]; then
                echo "Skipping decimation if nchan != 1 for all nsub. nsub=\${nsub##*_} nchan=\${nchan}"
                continue
            fi

            # Make full stokes and/or polarisation scrunched
            for stokes in ${meta.npols.join(' ')}; do
                if [ \${stokes} -eq 1 ]; then
                    # Polarisation scrunch option
                    stokes_op="-p"
                else
                    stokes_op=""
                fi

                echo "Decimate nsub=\${nsub##*_} nchan=\${nchan} stokes=\${stokes}"
                pam --setnsub \${nsub##*_} --setnchn \${nchan} -S \${stokes_op} -e \${nchan}ch_\${stokes}p_\${nsub}t.ar \${archive}
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
