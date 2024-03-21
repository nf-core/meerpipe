process PSRADD_CALIBRATE_CLEAN {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}", mode: params.publish_dir_mode, pattern: "${ template.baseName == "no_template" ? "*raw.ar" : "*zap.ar" }"

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(ephemeris), path(template)

    output:
    tuple val(meta), path(ephemeris), path(template), path("${meta.pulsar}_${meta.utc}_raw.ar"), path("${meta.pulsar}_${meta.utc}_zap.ar"), env(SNR), env(FLUX)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    """
    raw_only=${ template.baseName == "no_template" ? "true" : "false" }

    if ${params.use_edge_subints}; then
        # Grab all archives
        archives=\$(ls ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/*.ar)
    else
        if [ -z \$(ls ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/*.ar | head -n-1 | tail -n+2) ]; then
            # Grab all archives anyway because there are only two
            archives=\$(ls ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/*.ar)
        else
            # Grab all archives except for the first and last one
            archives=\$(ls ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/*.ar | head -n-1 | tail -n+2)
        fi
    fi

    if [ "\$raw_only" == "false" ]; then
        echo "Check if you need to change the template bins"
        obs_nbin=\$(vap -c nbin \$(echo "\$archives" | awk '{print \$1}') | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        std_nbin=\$(vap -c nbin ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        if [ "\$obs_nbin" == "\$std_nbin" ]; then
            std_template=${template}
        else
            echo "Making a new template with right number of bins"
            pam -b \$((std_nbin / obs_nbin)) -e new_std ${template}
            std_template=*new_std
        fi
    fi

    echo "Combine the archives"
    psradd \\
        -E ${ephemeris} \\
        -o ${meta.pulsar}_${meta.utc}_raw.ar \\
        \${archives}
    if [ "\$raw_only" == "false" ]; then
        echo "Clean the archive"
        clean_archive.py \\
            -a ${meta.pulsar}_${meta.utc}_raw.ar \\
            -T \${std_template} \\
            -o ${meta.pulsar}_${meta.utc}_zap.ar
    fi


    echo "Calibrate the polarisation of the archive"
    if [[ "${meta.cal_loc}" == "" || "${meta.cal_loc}" == "None" ]]; then
        # The archives have already be calibrated so just update the headers
        pac_args="-XP -e scalP"
    else
        # Use the Stokes paramaters files to calibrate the archive
        pac_args="-Q ${meta.cal_loc} -e scal"
    fi
    pac \${pac_args} -O ./ ${meta.pulsar}_${meta.utc}_raw.ar
    if [ "\$raw_only" == "false" ]; then
        pac \${pac_args} -O ./ ${meta.pulsar}_${meta.utc}_zap.ar
    fi

    echo "Update the RM value if available"
    rm_cat=\$(python -c "from meerpipe.data_load import RM_CAT;print(RM_CAT)")
    if grep -q "${meta.pulsar}" \${rm_cat}; then
        rm=\$(grep ${meta.pulsar} \${rm_cat} | tr -s ' ' | cut -d ' ' -f 2)
        echo "Found RM of \${rm} in the private RM catalogue"
    else
        rm=\$(psrcat -c RM ${meta.pulsar} -X -all | tr -s ' ' | cut -d ' ' -f 1)
        if [[ "\${rm}" == "*" || "\${rm}" == "WARNING:" ]]; then
            echo "No RM found in the ATNF catalogue"
            rm=0
        else
            echo "Found RM of \${rm} in the ATNF catalogue"
        fi
    fi
    pam --RM \${rm} -m ${meta.pulsar}_${meta.utc}_raw.scalP
    if [ "\$raw_only" == "false" ]; then
        pam --RM \${rm} -m ${meta.pulsar}_${meta.utc}_zap.scalP
    fi

    if \$raw_only; then
        echo "Delay correct"
        dlyfix -e ar ${meta.pulsar}_${meta.utc}_raw.scalP

        echo "Flux calibrate"
        # Create a time and polarisation scruchned profile
        pam -Tp -e tp ${meta.pulsar}_${meta.utc}_raw.ar
        fluxcal_meerkat \\
            --psr_name ${meta.pulsar} \\
            --obs_name ${meta.utc} \\
            --obs_header ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/obs.header \\
            --archive_file ${meta.pulsar}_${meta.utc}_raw.ar \\
            --tp_file *tp \\
            --par_file ${ephemeris}

        echo "No template provided so there will be no cleaned archive"
        touch ${meta.pulsar}_${meta.utc}_zap.ar
        SNR=None
        FLUX=None
    else
        echo "Delay correct"
        dlyfix -e ar ${meta.pulsar}_${meta.utc}_zap.scalP

        echo "Flux calibrate"
        # Create a time and polarisation scruchned profile
        pam -Tp -e tp ${meta.pulsar}_${meta.utc}_zap.ar
        fluxcal_meerkat \\
            --psr_name ${meta.pulsar} \\
            --obs_name ${meta.utc} \\
            --obs_header ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/obs.header \\
            --archive_file ${meta.pulsar}_${meta.utc}_zap.ar \\
            --tp_file *tp \\
            --par_file ${ephemeris}

        echo "Get the signal-to-noise ratio and flux density of the cleaned archive"
        pam -FTp -e FTp ${meta.pulsar}_${meta.utc}_zap.ar
        SNR=\$(psrstat -c snr=pdmp -c snr ${meta.pulsar}_${meta.utc}_zap.FTp | cut -d '=' -f 2)
        FLUX=\$(pdv -f ${meta.pulsar}_${meta.utc}_zap.FTp | tail -n 1 | tr -s ' ' | cut -d ' ' -f 7)
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${meta.pulsar}_${meta.utc}_raw.ar"
    touch "${meta.pulsar}_${meta.utc}_zap.ar"
    SNR=20
    """
}
