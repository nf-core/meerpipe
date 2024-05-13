process GRAB_ALL_PAIRS {
    tag "$pulsar"
    label 'process_single'
    label 'ephem_template'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:3.0.6' }"

    input:
    val pulsar

    output:
    path "all_pairs.csv", emit: out_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${pulsar}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    """
    if [[ -z "${params.ephemeris}" && -z "${params.template}" ]]; then
        grab_all_pairs ${pulsar} --out_dir ./
    else
        # Make directories
        mkdir -p ${params.project}/LBAND

        # Handling ephmeris
        if [ -z "${params.ephemeris}" ]; then
            echo "Grabbing ephemeris from repo"
            ephemeris="\$(grab_ephemeris ${pulsar} -p ${params.project})"
        else
            echo "Using input ephemeris"
            ephemeris="${params.ephemeris}"
        fi
        cp \${ephemeris} ${params.project}/${pulsar}.par
        ephemeris="${params.project}/${pulsar}.par"

        # Handling template
        if [ -z "${params.template}" ]; then
            echo "Grabbing template from repo"
            template="\$(grab_template ${pulsar} -p ${params.project})"
        else
            echo "Using input template"
            template="${params.template}"
        fi
        cp \${template} ${params.project}/LBAND/${pulsar}.std
        template="${params.project}/LBAND/${pulsar}.std"

        # Finally make output file
        echo "${pulsar},${params.project},\$(pwd)/\${ephemeris},\$(pwd)/\${template}" > all_pairs.csv
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch all_pairs.csv
    """
}
