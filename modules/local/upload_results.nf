process UPLOAD_RESULTS {
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
    tuple val(meta), path(png_files), path(dat_files), path(dynspec_files), path(results_json)

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
    import json
    import time
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, decode_id, get_graphql_id
    from psrdb.tables.pipeline_image import PipelineImage
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.toa import Toa

    logger = setup_logging(console=True, level=logging.DEBUG)
    if ${task.attempt} > 1:
        wait_time = (${task.attempt} - 1) * 30
        logger.info(f"Waiting for \${wait_time} s before retrying")
        time.sleep(wait_time)

    raw_only = ${ dynspec_files.baseName == "empty" ? "True" : "False" }

    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    pipeline_image_client = PipelineImage(client)
    toa_client            = Toa(client)
    pipeline_run_client   = PipelineRun(client)
    pipeline_run_client.get_dicts = True

    image_data = []
    # grab toa files
    for toa_file in glob("toa*png"):
        if "dmcorrected" in toa_file:
            type = "toa-dm-corrected"
        else:
            type = "toa-single"
        # file_loc, file_type, file_res, cleaned
        image_data.append( (toa_file, type, 'high', True) )

    # file_loc, file_type, file_res, cleaned
    image_data.append( ("raw_profile_ftp.png",    'profile',     'high', False) )
    image_data.append( ("raw_profile_fts.png",    'profile-pol', 'high', False) )
    image_data.append( ("raw_phase_time.png",     'phase-time',  'high', False) )
    image_data.append( ("raw_phase_freq.png",     'phase-freq',  'high', False) )
    image_data.append( ("raw_bandpass.png",       'bandpass',    'high', False) )
    image_data.append( ("raw_SNR_cumulative.png", 'snr-cumul',   'high', False) )
    image_data.append( ("raw_SNR_single.png",     'snr-single',  'high', False) )
    if os.path.exists("${meta.pulsar}_${meta.utc}_raw.ar.dynspec.png"):
        image_data.append( ("${meta.pulsar}_${meta.utc}_raw.ar.dynspec.png", 'dynamic-spectrum', 'high', False) )
    if not raw_only:
        image_data.append( ("cleaned_profile_ftp.png",    'profile',     'high', True ) )
        image_data.append( ("cleaned_profile_fts.png",    'profile-pol', 'high', True ) )
        image_data.append( ("cleaned_phase_time.png",     'phase-time',  'high', True ) )
        image_data.append( ("cleaned_phase_freq.png",     'phase-freq',  'high', True ) )
        image_data.append( ("cleaned_bandpass.png",       'bandpass',    'high', True ) )
        image_data.append( ("cleaned_SNR_cumulative.png", 'snr-cumul',   'high', True ) )
        if os.path.exists("${meta.pulsar}_${meta.utc}_zap.ar.dynspec.png"):
            image_data.append( ("${meta.pulsar}_${meta.utc}_zap.ar.dynspec.png", 'dynamic-spectrum', 'high', True ) )

    # Upload images
    for image_path, image_type, resolution, cleaned in image_data:
        image_response = pipeline_image_client.create(
            ${meta.pipe_id},
            image_path,
            image_type,
            resolution,
            cleaned,
        )
        content = json.loads(image_response.content)
        if image_response.status_code not in (200, 201):
            logger.error("Failed to upload image")
            exit(1)

    # Read in results JSON
    with open("results.json", "r") as f:
        results_dict = json.load(f)
    # Update pipeline run as completed (will update pulsarFoldResult)
    pipeline_run_response = pipeline_run_client.update(
        ${meta.pipe_id},
        "Completed",
        results_dict=results_dict,
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch done
    """
}
