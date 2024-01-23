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
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
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

    echo "Combine the archives and clean the archives in 200 file (~25 min) chunks"
    echo \${archives} | xargs -n 200 | while read -r archives_chunk; do
        first_file=\$(echo "\$archives_chunk" | awk '{print \$1}')
        last_file=\$(echo "\$archives_chunk" | awk '{print \$NF}')
        chunk_name="\$(basename \$first_file)_\$(basename \$last_file)"
        echo "Combine the archives of chunk \${chunk_name//.ar}"
        psradd \\
            -E ${ephemeris} \\
            -o ${meta.pulsar}_\${chunk_name//.ar}_chunk.raw \\
            \${archives_chunk}
        if [ "\$raw_only" == "false" ]; then
            echo "Clean the archive of chunk \${chunk_name//.ar}"
            ls -lrt
            clean_archive.py \\
                -a ${meta.pulsar}_\${chunk_name//.ar}_chunk.raw \\
                -T \${std_template} \\
                -o ${meta.pulsar}_\${chunk_name//.ar}_chunk.zap
        fi
    done
    echo "Combine all the file chunks"
    psradd -E ${ephemeris} -o ${meta.pulsar}_${meta.utc}_raw.ar *_chunk.raw
    if [ "\$raw_only" == "false" ]; then
        psradd -E ${ephemeris} -o ${meta.pulsar}_${meta.utc}_zap.ar *_chunk.zap
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
        fluxcal \\
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
        fluxcal \\
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
