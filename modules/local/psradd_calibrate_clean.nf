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
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(ephemeris), path(template)

    output:
    tuple val(meta), path(ephemeris), path(template), path("${meta.pulsar}_${meta.utc}_raw.ar"), path("${meta.pulsar}_${meta.utc}_zap.ar"), env(SNR)

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

    echo "Combine the archives"
    psradd -E ${ephemeris} -o ${meta.pulsar}_${meta.utc}_raw.raw \${archives}

    echo "Calibrate the polarisation of the archive"
    if [[ "${meta.cal_loc}" == "" || "${meta.cal_loc}" == "None" ]]; then
        # The archives have already be calibrated so just update the headers
        pac -XP -O ./ -e scalP ${meta.pulsar}_${meta.utc}_raw.raw
    else
        # Use the Stokes paramaters files to calibrate the archive
        pac -Q ${meta.cal_loc} -O ./ -e scal ${meta.pulsar}_${meta.utc}_raw.raw
    fi

    echo "Delay correct"
    dlyfix -e ar ${meta.pulsar}_${meta.utc}_raw.scalP

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
    pam --RM \${rm} -m ${meta.pulsar}_${meta.utc}_raw.ar

    if [ "${template.baseName}" == "no_template" ]; then
        echo "No template provided so not cleaning archive"
        touch ${meta.pulsar}_${meta.utc}_zap.ar
        SNR=None
    else
        echo "Check if you need to change the template bins"
        obs_nbin=\$(vap -c nbin ${meta.pulsar}_${meta.utc}_raw.raw | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        std_nbin=\$(vap -c nbin ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        if [ "\$obs_nbin" == "\$std_nbin"]; then
            std_template=${template}
        else
            echo "Making a new template with right number of bins"
            pam -b \$((std_nbin / obs_nbin)) -e new_std ${template}
            std_template=*new_std
        fi
        echo "Clean the archive"
        clean_archive.py -a ${meta.pulsar}_${meta.utc}_raw.ar -T \${std_template} -o ${meta.pulsar}_${meta.utc}_zap.ar

        # Get the signal to noise ratio of the cleaned archive
        SNR=\$(psrstat -j FTp -c snr=pdmp -c snr ${meta.pulsar}_${meta.utc}_zap.ar | cut -d '=' -f 2)

        echo "Flux calibrate"
        # Create a time and polarisation scruchned profile
        pam -Tp -e tp ${meta.pulsar}_${meta.utc}_zap.ar
        fluxcal -psrname ${meta.pulsar} -obsname ${meta.utc} -obsheader ${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/obs.header -cleanedfile ${meta.pulsar}_${meta.utc}_zap.ar -rawfile ${meta.pulsar}_${meta.utc}_raw.ar -tpfile *tp -parfile ${ephemeris}
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
