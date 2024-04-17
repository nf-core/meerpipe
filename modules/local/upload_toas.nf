process UPLOAD_TOAS {
    tag "$meta.id"
    label 'psrdb'
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:latest' }"

    input:
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
    template_band = template.split("/")[-2]
    template_project = template.split("/")[-3]
    template_response = template_client.create(
        "${meta.pulsar}",
        template_band,
        template,
        project_short=template_project,
    )
    logger.debug(template_response)
    template_id = get_rest_api_id(template_response, logging.getLogger(__name__))

    # Upload TOAs
    for nchan in [${meta.nchans.join(', ')}]:
        for npol in [${meta.npols.join(', ')}]:
            # Grab all the nsub types
            nsub_types = []
            if "${params.use_all_nsub}" == "true":
                nsub_types.append("all")
            if "${params.use_mode_nsub}" == "true":
                nsub_types.append("mode")
            if "${params.use_max_nsub}" == "true":
                nsub_types.append("max")

            # Work out which types only have 1 nsub
            nsub_1_types = []
            nsub_other_types = []
            for nsub_type in nsub_types:
                print(f"${meta.pulsar}_*{nchan}ch_{npol}p_{nsub_type}*t.ar.tim")
                toa_files = glob(f"${meta.pulsar}_*{nchan}ch_{npol}p_{nsub_type}*t.ar.tim")
                if len(toa_files) == 1:
                    toa_file = toa_files[0]
                    nsub = int(toa_file.split(f"{npol}p_{nsub_type}_")[-1].split("t.ar.tim")[0])
                    if nsub == 1:
                        nsub_1_types.append(nsub_type)
                    else:
                        nsub_other_types.append((nsub_type, toa_file))

            # Upload the nsub=1 type using just the nsub=1 file
            toa_file = glob(f"${meta.pulsar}_*{nchan}ch_{npol}p_1t.ar.tim")[0]
            all_nsubs = True if "all" in nsub_1_types else False
            max_nsubs = True if "max" in nsub_1_types else False
            mode_nsubs = True if "mode" in nsub_1_types else False
            logger.info(f"Uploading Toa nsub=1 file {toa_file} with nchan={nchan}, npol={npol}, all_nsubs={all_nsubs}, max_nsubs={max_nsubs}, mode_nsubs={mode_nsubs},")
            with open(toa_file, "r") as f:
                toa_lines = f.readlines()
                toa_response = toa_client.create(
                    ${meta.pipe_id},
                    "${meta.project_short}",
                    "${ephemeris}",
                    template_id,
                    toa_lines,
                    dmCorrected=False,
                    minimumNsubs=True,
                    maximumNsubs=max_nsubs,
                    allNsubs=all_nsubs,
                    modeNsubs=mode_nsubs,
                    npol=npol,
                    nchan=nchan,
                )
                if toa_response.status_code not in (200, 201):
                    logger.error("Failed to upload TOA")
                    exit(1)
                logger.info(get_graphql_id(toa_response, "toa", logger))

            # Upload all the non nsub=1 types
            for nsub_type, toa_file in nsub_other_types:
                all_nsubs = True if "all" == nsub_type else False
                max_nsubs = True if "max" == nsub_type else False
                mode_nsubs = True if "mode" == nsub_type else False
                logger.info(f"Uploading Toa file {toa_file} with nchan={nchan}, npol={npol}, all_nsubs={all_nsubs}, max_nsubs={max_nsubs}, mode_nsubs={mode_nsubs},")
                with open(toa_file, "r") as f:
                    toa_lines = f.readlines()
                    toa_response = toa_client.create(
                        ${meta.pipe_id},
                        "${meta.project_short}",
                        "${ephemeris}",
                        template_id,
                        toa_lines,
                        dmCorrected=False,
                        minimumNsubs=False,
                        maximumNsubs=max_nsubs,
                        allNsubs=all_nsubs,
                        modeNsubs=mode_nsubs,
                        npol=npol,
                        nchan=nchan,
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
