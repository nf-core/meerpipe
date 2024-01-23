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

process UPLOAD_TOAS {
    tag "$meta.id"
    label 'psrdb'
    label 'process_single'

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
    tuple val(meta), path(ephemeris), path(template), path(toas)

    output:
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
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #!/usr/bin/env python

    import os
    import time
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_graphql_id, get_rest_api_id
    from psrdb.tables.toa import Toa
    from psrdb.tables.template import Template

    logger = setup_logging(console=True, level=logging.DEBUG)
    if ${task.attempt} > 1:
        wait_time = (${task.attempt} - 1) * 30
        logger.info(f"Waiting for \${wait_time} s before retrying")
        time.sleep(wait_time)

    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    toa_client      = Toa(client)
    template_client = Template(client)

    # Upload template (if not uploaded already)
    template = os.path.realpath("${template}")
    template_band = template.split("/")[-3]
    template_project = template.split("/")[-2]
    template_response = template_client.create(
        "${meta.pulsar}",
        template_band,
        template,
        project_short=template_project,
    )
    logger.debug(template_response)
    template_id = get_rest_api_id(template_response, logging.getLogger(__name__))

    # Upload TOAs
    for toa_file in ["${toas.join('","')}"]:
        if "dm_corrected" in toa_file:
            dmcorrected = True
            # TODO don't skip this if it's useful
            continue
        else:
            dmcorrected = False

        # Work out if this is the minimum or maximum number of subints
        nchan = toa_file.split("_zap.")[-1].split("ch")[0]
        nsub = toa_file.split("_zap."+nchan+"ch1p")[-1].split("t.ar")[0]
        toas_same_nchan = glob("*_zap." + nchan + "ch*" + toa_file.split(".ar.")[1])
        nsubs_list = []
        for toa in toas_same_nchan:
            nsubs_list.append(int(toa.split("_zap."+nchan+"ch1p")[-1].split("t.ar")[0]))
        minimum_nsubs = False
        maximum_nsubs = False
        if max(nsubs_list) == int(nsub):
            maximum_nsubs = True
        if min(nsubs_list) == int(nsub):
            minimum_nsubs = True

        logger.info(f"Uploading Toa file {toa_file} with maximum_nsubs={maximum_nsubs} and minimum_nsubs={minimum_nsubs}")
        with open(toa_file, "r") as f:
            toa_lines = f.readlines()
            toa_response = toa_client.create(
                ${meta.pipe_id},
                "${meta.project_short}",
                "${ephemeris}",
                template_id,
                toa_lines,
                dmcorrected,
                minimum_nsubs,
                maximum_nsubs,
            )
            if toa_response.status_code not in (200, 201):
                logger.error("Failed to upload TOA")
                exit(1)
            logger.info(get_graphql_id(toa_response, "toa", logger))
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch done
    """
}
