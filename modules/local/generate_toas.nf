process GENERATE_TOAS {
    tag "$meta.id"
    label 'process_medium'
    label 'psrchive'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/timing/${meta.project_short}", mode: 'copy', pattern: "*.{tim,par,std}"


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(ephemeris), path(template), path(decimated_archives)

    output:
    tuple val(meta), path(ephemeris), path(template), path("*.tim")

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
    # Loop over each DECIMATEd archive
    for ar in ${decimated_archives.join(' ')}; do
        # Grab archive and template nchan
        nchan=\$(vap -c nchan \$ar | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        tnchan=\$(vap -c nchan ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)

        # Use portrait mode if template has more frequency channels
        if [ "\$tnchan" -gt "\$nchan" ]; then
            port="-P"
        else
            port=""
        fi

        echo "Generating TOAs for \${ar}.tim\n----------------------------------"
        pat -jp \$port  -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM \$ar  > \${ar}.tim
    done
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${meta.pulsar}_${meta.utc}_zap.1ch1p1t.tim
    """
}
