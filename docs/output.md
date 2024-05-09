# nf-core/meerpipe: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Meerpipe products](#meerpipe-products) - Raw read QC
  - [Images](#images) - Generate images for the MeerTime data portal
  - [Decimated](#decimated) - Decimate the data for timing analysis
  - [Timing](#timing) - Create ToA files for timing analysis
  - [Scintillation](#scintillation) - Create files for scintillation analysis
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Meerpipe products

The processed fold mode data is generated with `meerpipe` using the raw fold mode data and is stored in the following directory:

```
${params.outdir}/<pulsar>/<utc>/<beam>
```

This directory will contain a `results.json` file which are outputs of `meerpipe` calculations
and either a cleaned/zapped (`<pulsar>_<utc>_zap.ar`) or raw (`<pulsar>_<utc>_raw.ar`) file if no template is available (templates are required for cleaning with [MeerGuard](https://github.com/danielreardon/MeerGuard)).
If there is only a raw file you can either manually clean it with
[paz](https://ozgrav.github.io/2023-09-25_NWU_Pulsar_Timing_Workshop/PulsarData/index.html#pazi)
or [pazi](https://ozgrav.github.io/2023-09-25_NWU_Pulsar_Timing_Workshop/PulsarData/index.html#paz)
or add a template to the [ephemeris and template repo](/meerkat_pulsar_docs/ephem_template/#development-add-or-update-ephemerides-and-templates) so that the `meerpipe` pipeline can clean it and create the following outputs.

This directory also contains the following subdirectories that will be explained in the following subsections:

- [`images`](#images)
- [`decimated`](#decimated)
- [`timing`](#timing)
- [`scintillation`](#scintillation)

### Images

The images directory (`${params.outdir}/<pulsar>/<utc>/<beam>/images/`) contains all the images that will be uploaded to the [MeerTime data portal](https://pulsars.org.au/).
Each image will either start with `cleaned` or `raw` depending on whether the image was created from a cleaned/zapped or raw archive respectively.
The images names and descriptions are listed below in the same order they appear on the data portal:

- `{cleaned/raw}_profile_ftp.png`: The polarisation pulse profile.
- `{cleaned/raw}_profile_fts.png`: The pulse profile.
- `{cleaned/raw}_phase_time.png`: The phase vs. time plot.
- `{cleaned/raw}_phase_freq.png`: The phase vs. frequency plot.
- `{cleaned/raw}_bandpass.png`: The bandpass plot which shows the frequency response and which channels were flagged.
- `{cleaned/raw}_SNR_cumulative.png`: The cumulative signal-to-noise ratio plot which shows how the SNR increases with time.
- `{cleaned/raw}_SNR_single.png`: The single subint signal-to-noise ratio plot which shows the SNR at each subint.
- `cleaned_rmfit.png`: The result of the `rmfit` command for checking the quality of the RM measurement.

### Decimated

The decimated directory (`${params.outdir}/<pulsar>/<utc>/<beam>/decimated/`) contains the decimated archives that are used for the timing analysis for all projects.
The decimated archives have the following naming convention:

```
<pulsar>_<utc>_{zap/raw}<_chopped>.<nchan>ch_<npol>p_<nsub_type>_<nsub>t.ar
```

where `{zap/raw}` represents if the archive was cleaned or not,
`<_chopped>` is included in the file name if the edge frequency channels are removed,
`<nchan>` is the number of frequency channels,
`<npol>` is the number of polarisations,
`<nsub_type>` is the method used to calculate the number of subintegrations (see [next section](#nsub-types)), and
`<nsub>` is the number of time subintegrations.

The following are examples of the decimated file names:

- `J1744-1134_2019-10-05-11:17:35_zap_chopped.16ch_1p_1t.ar` is a decimated cleaned archive with the edge channels removed, 16 frequency channels, 1 polarisation and 1 subintegration. Note that the `1` nsub type observations don't have both their `<nsub_type>` and `<nsub>` as they are both the same.
- `J1744-1134_2019-10-05-11:17:35_zap_chopped.1ch_4p_all_32t.ar` is a decimated cleaned archive with the edge channels removed, 1 frequency channel (frequency scrunched), full Stokes polarisations and all subintegrations (not time scrunched) which is 32 for this observation.
- `J1744-1134_2019-10-05-11:17:35_raw_chopped.928ch_4p_mode_1t.ar` is a decimated RAW archive with the edge channels removed, 928 frequency channels, full Stokes polarisations and mode subintegrations which is 1 for this observation. Note this is a raw archive so caution should be used when using it for sensitive science.

#### Nsub types

There are currently four methods of how to calculate how many time subintegrations to use for an observations which are listed below.

- "1": a single nsub (time scrunched)
- "all": all available nsubs (no time scrunching), only done for single nchan decimations (frequency scrunched)
- "max" the maximum number of subints possible for the observation based on the S/N ratio.
The maximum is calculated using the `meerpipe` script [`calc_max_nsub`](https://github.com/OZGrav/meerpipe/blob/main/meerpipe/scripts/calc_max_nsub.py) which uses the logic:
$$
nsub = \left ( \frac{S/N}{S/N_D} \right )^2 \frac{1}{nchan}
$$
where $S/N$ is the signal-to-noise ration of the observation,
$S/N_D$ is the desired signal-to-noise ratio of ToA (12 by default) and
$nchan$ is the number of frequency channels for this decimation.
The minimum length of each of these sub-integrations must be 480 seconds (to prevent a huge number of ToAs).
- "mode" the length of each subintegration is equal to the most common observation duration.
This value is calculated as part of the [webportal query](https://gitlab.com/CAS-eResearch/GWDC/meertime_dataportal/-/blob/main/backend/dataportal/graphql/queries.py?ref_type=heads#L338) which rounds all values to the nearest 32 seconds and finds the most common duration, prioritising short observations in the case of a draw.

### Timing

The timing directory contains subdirectories for each project (e.g PTA) (`${params.outdir}/<pulsar>/<utc>/<beam>/decimated/<project>/`).
In each of these project timing subdirectories there is a ephemeris file (`<pulsar>.par`) and a template file (`<pulsar>.std`) which are used to make a time of arrival (ToA) `.tim` files.
The ToA files have the following format:

```
<pulsar>_<utc>_zap<_chopped>.<nchan>ch_<npol>p_<nsub_type>_<nsub>.tim
```

where, as above, `<_chopped>` is included in the file name if the edge frequency channels are removed,
`<nchan>` is the number of frequency channels,
`<npol>` is the number of polarisations,
`<nsub_type>` is the method used to calculate the number of subintegrations (see [previous section](#nsub-types)), and
`<nsub>` is the number of time subintegrations.

The ToAs can manually be combined into a single `.tim` file or more easily downloaded using `psrdb toa download`, see the [psrd docs](https://psrdb.readthedocs.io/en/latest/how_to_use.html#toa-download-example) for examples of how to do so.



### Scintillation

The scintillation directory (`${params.outdir}/<pulsar>/<utc>/<beam>/scintillation/`) contains files used for scintillation analysis.
The `.dynspec` data files are generated with the `psrflux` psrchive script and have the following naming convention:

```
<pulsar>_<utc>_{raw,zap}.ar.dynspec
```

Where `raw` is for the raw archives and `zap` is for the cleaned/zapped archives.

There are also png files which are created with the [scintools](https://github.com/danielreardon/scintools) repository and have the following naming convention:

```
<pulsar>_<utc>_{raw,zap}.ar.dynspec.png
```

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
