process OBS_LIST {
    label 'process_single'
    label 'psrdb'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:latest' }"

    input:
    val pulsars
    val utcs
    val utce
    val project_short
    path obs_csv
    val upload
    val psrdb_url
    val psrdb_token
    val ephemeris
    val template
    val outdir

    output:
    path "processing_jobs.csv", emit: out_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import os
    import re
    import json
    import time
    import base64
    import shutil
    import logging
    import pandas as pd
    from datetime import datetime

    from psrdb.tables.observation import Observation
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.template import Template
    from psrdb.tables.ephemeris import Ephemeris
    from psrdb.tables.pulsar_fold_result import PulsarFoldResult
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_rest_api_id, get_graphql_id, decode_id
    from ephem_template.python_grabber import grab_ephemeris, grab_template, NoFilesFound

    # Get pipeline config from the params and manifest
    params_string = '${params.all()}'
    pipeline_config = {}
    # Extracting key-value pairs from the manifest
    manifest_match = re.search(r'--manifest "(.*?)"', params_string)
    if manifest_match:
        manifest_string = manifest_match.group(1)
        manifest_list = re.findall(r'(\\w+):([^,]+)', manifest_string)
        manifest_dict = dict(manifest_list)
        pipeline_config['manifest'] = manifest_dict
        # Remove the manifest string from the params string
        params_string = params_string.replace(manifest_string, "").replace("--manifest", "")
    # Extracting key-value pairs from the params string
    pairs = params_string.split('--')
    for pair in pairs:
        if pair.strip():  # Skip empty strings
            key_value = pair.strip().split(maxsplit=1)
            if len(key_value) == 2:
                key, value = key_value
            else:
                key = key_value[0]
                value = ""
            pipeline_config[key.strip()] = value.strip()

    # PSRDB setup
    logger = setup_logging(level=logging.DEBUG)
    client = GraphQLClient("${psrdb_url}", "${psrdb_token}", logger=logger)
    pipe_run_client  = PipelineRun(client)
    template_client  = Template(client)
    ephemeris_client = Ephemeris(client)
    obs_client       = Observation(client)
    obs_client.get_dicts = True
    obs_client.set_use_pagination(True)
    pfr_client       = PulsarFoldResult(client)
    pfr_client.get_dicts = True
    pfr_client.set_use_pagination(True)

    if "${obs_csv.baseName}" == "none_given":
        if "${pulsars}" in ("", "null"):
            pulsar_list = ""
        else:
            pulsar_list = ["${pulsars.split(',').join('","')}"]
        print(pulsar_list)
        # Query based on provided parameters
        obs_data = obs_client.list(
            pulsar_name=pulsar_list,
            project_short="${project_short}",
            main_project="MeerTIME",
            utcs="${utcs}",
            utce="${utce}",
            obs_type="fold",
        )
        obs_df = pd.DataFrame(columns=["Obs ID","Pulsar Jname","UTC Start","Project Short Name","Beam #","Observing Band","Duration (s)","Mode Duration (s)","Nchan","Nbin","Calibration Location"])
        for obs in obs_data:
            # if obs['pulsar']['name'] in ["J1735-3028A"]:
            #     continue
            obs_df = pd.concat(
                [
                    obs_df,
                    pd.Series({
                        "Obs ID": decode_id(obs['id']),
                        "Pulsar Jname": obs['pulsar']['name'],
                        "UTC Start": datetime.strptime(
                            obs['utcStart'],
                            '%Y-%m-%dT%H:%M:%S+00:00',
                        ).strftime('%Y-%m-%d-%H:%M:%S'),
                        "Project Short Name": obs['project']['short'],
                        "Beam #": obs['beam'],
                        "Observing Band": obs['band'],
                        "Duration (s)": obs['duration'],
                        "Mode Duration (s)": obs['modeDuration'],
                        "Nchan": obs["foldNchan"],
                        "Nbin": obs["foldNbin"],
                        "Calibration Location": obs['calibration']['location'],
                    }).to_frame().T,
                ],
                ignore_index=True
            )
    else:
        # Read in obs from csv
        obs_df = pd.read_csv("${obs_csv}")


    # Grab all the ephems and templates for each pulsar first
    pulsar_ephem_template = {}
    for _, obs in obs_df.iterrows():
        # Extract data from obs_df
        pulsar   = obs['Pulsar Jname']
        project  = obs['Project Short Name']
        band     = obs['Observing Band']
        if pulsar in pulsar_ephem_template.keys():
            if band in pulsar_ephem_template[pulsar].keys():
                # Already grabbed so continue
                continue
        else:
            pulsar_ephem_template[pulsar] = {}

        # Grab ephermis and templates
        if "${ephemeris}" in ("", "null"):
            ephemeris = grab_ephemeris(pulsar, project, fold=True)
            os.makedirs(project, exist_ok=True)
            shutil.copyfile(ephemeris, f"{project}/{pulsar}.par")
            ephemeris = os.path.join(os.getcwd(), f"{project}/{pulsar}.par")
        else:
            os.makedirs("${params.project}", exist_ok=True)
            shutil.copyfile("${ephemeris}", f"${params.project}/{pulsar}.par")
            ephemeris = os.path.join(os.getcwd(), f"${params.project}/{pulsar}.par")
        if "${template}" in ("", "null"):
            try:
                template = grab_template(pulsar, project, band, fold=True)
                os.makedirs(f"{band}", exist_ok=True)
                os.makedirs(f"{band}/{project}", exist_ok=True)
                shutil.copyfile(template, f"{band}/{project}/{pulsar}.std")
                template = os.path.join(os.getcwd(), f"{band}/{project}/{pulsar}.std")
            except NoFilesFound:
                template = os.path.join(os.getcwd(), "no_template.std")
                with open(template, 'w'):
                    pass # Make an empty file
        else:
            os.makedirs("LBAND", exist_ok=True)
            os.makedirs("LBAND/${params.project}", exist_ok=True)
            shutil.copyfile("${template}", f"LBAND/${params.project}/{pulsar}.std")
            template = os.path.join(os.getcwd(), f"LBAND/${params.project}/{pulsar}.std")
        ephem_template = {
            "ephemeris": ephemeris,
            "template": template,
        }
        if "${upload}" == "true":
            # Get or create template
            if template == os.path.join(os.getcwd(), "no_template.std"):
                ephem_template["template_id"] = -1
            else:
                template_band = template.replace("/fold/", "/").split("/")[-3]
                template_project = template.replace("/fold/", "/").split("/")[-2]
                if template_project == "unknown":
                    template_project = "PTUSE"
                logger.info(f"Creating template with band: {template_band} project: {template_project} template: {template}")
                template_response = template_client.create(
                    pulsar,
                    template_band,
                    template,
                    project_short=template_project,
                )
                logger.debug(template_response)
                ephem_template["template_id"] = get_rest_api_id(template_response, logging.getLogger(__name__))
            # Get or create ephemeris
            ephemeris_project = ephemeris.replace("/fold/", "/").split("/")[-2]
            if ephemeris_project == "unknown":
                ephemeris_project = "PTUSE"
            logger.info(f"Creating ephemeris with project: {ephemeris_project} ephemeris: {ephemeris}")
            ephemeris_response = ephemeris_client.create(
                pulsar,
                ephemeris,
                project_short=ephemeris_project,
                comment="",
            )
            logger.debug(ephemeris_response)
            ephem_template["ephemeris_id"] = get_graphql_id(ephemeris_response, "ephemeris", logging.getLogger(__name__))
        pulsar_ephem_template[pulsar][band] = ephem_template



    # Add ephemeris and template to df and create the pipeline run object
    obs_df['pipe_id'] = ''
    obs_df['ephemeris'] = ''
    obs_df['template'] = ''
    obs_df['n_obs'] = 0
    for index, obs in obs_df.iterrows():
        pulsar = obs['Pulsar Jname']
        band   = obs['Observing Band']
        ephemeris = pulsar_ephem_template[pulsar][band]["ephemeris"]
        template  = pulsar_ephem_template[pulsar][band]["template"]

        logger.info(f"Setting up {index+1}/{len(obs_df)} ID: {obs['Obs ID']} pulsar: {pulsar} band: {band} template: {template} ephemeris: {ephemeris}")

        # Set job as running
        if "${upload}" == "true":
            ephemeris_id = pulsar_ephem_template[pulsar][band]["ephemeris_id"]
            template_id  = pulsar_ephem_template[pulsar][band]["template_id"]

            try:
                pipe_run_data = pipe_run_client.create(
                    obs['Obs ID'],
                    ephemeris_id,
                    template_id,
                    pipeline_config['manifest']["name"],
                    pipeline_config['manifest']["description"],
                    pipeline_config['manifest']["version"],
                    "Running",
                    "${outdir}",
                    pipeline_config,
                )
                pipe_id = get_graphql_id(pipe_run_data, "pipeline_run", logging.getLogger(__name__))
            except ValueError as e:
                logger.error(f"Failed to create pipeline run for {obs['Obs ID']}: {e}")
                logger.info(f"Waiting 30 seconds and trying again")
                time.sleep(30)
                pipe_run_data = pipe_run_client.create(
                    obs['Obs ID'],
                    ephemeris_id,
                    template_id,
                    pipeline_config['manifest']["name"],
                    pipeline_config['manifest']["description"],
                    pipeline_config['manifest']["version"],
                    "Running",
                    "${outdir}",
                    pipeline_config,
                )
                pipe_id = get_graphql_id(pipe_run_data, "pipeline_run", logging.getLogger(__name__))
        else:
            # No uploading so don't make a processing item
            pipe_id = None
        obs_df.at[index, 'pipe_id']   = pipe_id
        obs_df.at[index, 'ephemeris'] = ephemeris
        obs_df.at[index, 'template']  = template

        if obs["Calibration Location"] in ("", "null", "None", None):
            obs_df.at[index, 'Calibration Location'] = os.path.join(os.getcwd(), "no_cal_file")
            with open("no_cal_file", 'w'):
                pass # Make an empty file

        # Count observations with the same pulsar
        obs_df.at[index, 'n_obs']     = len(obs_df[obs_df['Pulsar Jname'] == pulsar])

    if "${params.use_prev_ar}" == "true":
        obs_df['sn'] = 0.
        obs_df['flux'] = 0.
        obs_df['percent_rfi_zapped'] = 0.
        with open('empty_raw.ar', 'w'):
            pass # Make an empty file
        obs_df['raw_archive'] = os.path.join(os.getcwd(), 'empty_raw.ar')
        obs_df['clean_archive'] = ''
        for index, obs in obs_df.iterrows():
            pfr_data = pfr_client.list(
                pulsar=obs["Pulsar Jname"],
                mainProject="MeerTIME",
                utcStart=obs["UTC Start"],
                beam=obs["Beam #"],
            )
            print(pfr_data)
            obs_df.at[index, 'sn'] = pfr_data[0]['pipelineRun']['sn']
            obs_df.at[index, 'flux'] = pfr_data[0]['pipelineRun']['flux']
            obs_df.at[index, 'percent_rfi_zapped'] = pfr_data[0]['pipelineRun']['percentRfiZapped']
            obs_df.at[index, 'clean_archive'] = f"${params.outdir}/{obs['Pulsar Jname']}/{obs['UTC Start']}/{obs['Beam #']}/{obs['Pulsar Jname']}_{obs['UTC Start']}_zap.ar"

    # Write out results
    obs_df.drop('Obs ID', axis=1, inplace=True)
    obs_df.to_csv("processing_jobs.csv", header=False, index=False)
    """

    stub:
    """
    touch processing_jobs.csv
    """
}
