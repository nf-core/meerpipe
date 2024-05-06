process GENERATE_RESIDUALS {
    tag "$meta.id"
    label 'process_single'
    label 'psrdb_tempo2'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:latest' }"

    input:
    tuple val(meta), path(ephemeris)

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
    # Loop over each of the TOA filters set in params
    nsub_types="1"
    if ${params.use_max_nsub}; then
        nsub_types="\${nsub_types} max"
    fi
    if ${params.use_mode_nsub}; then
        nsub_types="\${nsub_types} mode"
    fi
    for nsub_type in \${nsub_types}; do
        for nchan in ${meta.nchans.join(' ')}; do
            if ((nchan > ${params.max_nchan_upload})); then
                echo -e "\\nSkipping the calculation and upload for residuals for \${nsub_type} \${nchan}\\n--------------------------\\n"
                continue
            fi
            echo -e "\\nDownload the toa file and fit the residuals for nsub=\${nsub_type} nchan=\${nchan} npol=1\\n--------------------------\\n"
            psrdb -u ${params.psrdb_url} -t ${params.psrdb_token} toa download ${meta.pulsar} --project ${meta.project_short} --nsub_type \$nsub_type --nchan \$nchan --npol 1
            echo -e "\\nGenerating residuals for nsub=\${nsub_type} nchan=\${nchan} npol=1\\n--------------------------\\n"
            tempo2_wrapper.sh toa_${meta.pulsar}_${meta.project_short}_\${nsub_type}_nsub_nchan\${nchan}_npol1.tim ${ephemeris}
            if [ -f "toa_${meta.pulsar}_${meta.project_short}_\${nsub_type}_nsub_nchan\${nchan}_npol1.tim.residual" ]; then
                echo -e "\\nUpload the residuals\\n--------------------------\\n"
                psrdb -u ${params.psrdb_url} -t ${params.psrdb_token} residual create toa_${meta.pulsar}_${meta.project_short}_\${nsub_type}_nsub_nchan\${nchan}_npol1.tim.residual
            else
                echo -e "\\nResiduals file toa_${meta.pulsar}_${meta.project_short}_\${nsub_type}_nsub_nchan\${nchan}_npol1.tim.residual not found\\n--------------------------\\n"
            fi
        done
    done
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch done
    """
}
