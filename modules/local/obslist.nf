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

process OBS_LIST {
    label 'psrdb'

    input:
    val pulsars
    val utcs
    val utce
    val project_short
    val obs_csv
    val upload
    val psrdb_url
    val psrdb_token
    val ephemeris
    val template
    val outdir
    path manifest

    output:
    path "processing_jobs.csv", emit: out_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import os
    import json
    import base64
    import logging
    import pandas as pd
    from datetime import datetime

    from psrdb.tables.observation import Observation
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.template import Template
    from psrdb.tables.ephemeris import Ephemeris
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_rest_api_id, get_graphql_id, decode_id
    from ephem_template.python_grabber import grab_ephemeris, grab_template

    # PSRDB setup
    logger = setup_logging(level=logging.DEBUG)
    client = GraphQLClient("${psrdb_url}", "${psrdb_token}", logger=logger)
    obs_client       = Observation(client)
    pipe_run_client  = PipelineRun(client)
    template_client  = Template(client)
    ephemeris_client = Ephemeris(client)
    obs_client.get_dicts = True
    obs_client.set_use_pagination(True)

    if "${obs_csv}" == "":
        if "${pulsars}" == "":
            pulsar_list = ""
        else:
            pulsar_list = ["${pulsars.split(',').join('","')}"]
        print(pulsar_list)
        # Query based on provided parameters
        obs_data = obs_client.list(
            pulsar_name=pulsar_list,
            project_short="${project_short}",
            utcs="${utcs}",
            utce="${utce}",
            obs_type="fold",
        )
        obs_df = pd.DataFrame(columns=["Obs ID","Pulsar Jname","UTC Start","Project Short Name","Beam #","Observing Band","Duration (s)","Calibration Location"])
        for obs in obs_data:
            if obs['pulsar']['name'] in ["J1735-3028A"]:
                continue
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
        if "${ephemeris}" == "":
            ephemeris = grab_ephemeris(pulsar, project)
        else:
            ephemeris = "${ephemeris}"
        if "${template}" == "":
            try:
                template = grab_template(pulsar, project, band)
            except ValueError:
                template = os.path.join(os.getcwd(), "no_template.std")
                with open(template, 'w'):
                    pass # Make an empty file
        else:
            template = "${template}"
        ephem_template = {
            "ephemeris": ephemeris,
            "template": template,
        }
        if "${upload}" == "true":
            # Get or create template
            if template == os.path.join(os.getcwd(), "no_template.std"):
                ephem_template["template_id"] = -1
            else:
                template_band = template.split("/")[-3]
                template_project = template.split("/")[-2]
                template_response = template_client.create(
                    pulsar,
                    template_band,
                    template,
                    project_short=template_project,
                )
                logger.debug(template_response)
                ephem_template["template_id"] = get_rest_api_id(template_response, logging.getLogger(__name__))
            # Get or create ephemeris
            ephemeris_project = ephemeris.split("/")[-2]
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
    for index, obs in obs_df.iterrows():
        pulsar = obs['Pulsar Jname']
        band   = obs['Observing Band']
        ephemeris = pulsar_ephem_template[pulsar][band]["ephemeris"]
        template  = pulsar_ephem_template[pulsar][band]["template"]

        logger.info(f"Setting up ID: {obs['Obs ID']} pulsar: {pulsar} band: {band} template: {template} ephemeris: {ephemeris}")

        # Set job as running
        if "${upload}" == "true":
            ephemeris_id = pulsar_ephem_template[pulsar][band]["ephemeris_id"]
            template_id  = pulsar_ephem_template[pulsar][band]["template_id"]
            with open("${manifest}", 'r') as file:
                # Load the JSON data
                pipeline_config = json.load(file)

            pipe_run_data = pipe_run_client.create(
                obs['Obs ID'],
                ephemeris_id,
                template_id,
                pipeline_config["pipeline_name"],
                pipeline_config["pipeline_description"],
                pipeline_config["pipeline_version"],
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

    # Write out results
    obs_df.drop('Obs ID', axis=1, inplace=True)
    obs_df.to_csv("processing_jobs.csv", header=False, index=False)
    """

    stub:
    """
    touch processing_jobs.csv
    """
}
