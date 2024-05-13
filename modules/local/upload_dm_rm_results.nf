process UPLOAD_DM_RM_RESULTS {
    tag "$meta.id"
    label 'psrdb'
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:3.0.6' }"

    input:
    tuple val(meta), path(results_json), path(png_files)

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

    def return_none_or_float(value):
        if value == "None":
            return None
        else:
            return float(value)

    logger = setup_logging(console=True, level=logging.DEBUG)
    if ${task.attempt} > 1:
        wait_time = (${task.attempt} - 1) * 30
        logger.info(f"Waiting for \${wait_time} s before retrying")
        time.sleep(wait_time)

    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    pipeline_image_client = PipelineImage(client)
    toa_client            = Toa(client)
    pipeline_run_client   = PipelineRun(client)
    pipeline_run_client.get_dicts = True

    image_data = []
    if os.path.exists("cleaned_rmfit.png"):
        image_data.append( ("cleaned_rmfit.png", 'rmfit', 'high', True ) )

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
    with open("${meta.pulsar}_${meta.utc}_dm_rm_fit.json", "r") as f:
        results_dict = json.load(f)
    # Add sn and flux from metadata
    results_dict["dm"]       = return_none_or_float(results_dict["DM"])
    results_dict["dm_err"]   = return_none_or_float(results_dict["ERR"])
    results_dict["dm_epoch"] = return_none_or_float(results_dict["EPOCH"])
    results_dict["dm_chi2r"] = return_none_or_float(results_dict["CHI2R"])
    results_dict["dm_tres"]  = return_none_or_float(results_dict["TRES"])
    results_dict["rm"]       = return_none_or_float(results_dict["RM"])
    results_dict["rm_err"]   = return_none_or_float(results_dict["RM_ERR"])
    results_dict["sn"]       = float(${meta.snr})
    results_dict["flux"]     = float(${meta.flux})
    results_dict["percent_rfi_zapped"] = float(${meta.percent_rfi_zapped})

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
