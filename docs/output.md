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


### Decimated

The decimated directory (`${params.outdir}/<pulsar>/<utc>/<beam>/decimated/`) contains the decimated archives that are used for the timing analysis for all projects.
The decimated archives have the following naming convention:

```
<pulsar>_<utc>_zap<_chopped>.<nchan>ch<npol>p<nsub>t.ar
```

where `<_chopped>` is included in the file name if the edge frequency channels are removed,
`<nchan>` is the number of frequency channels,
`<npol>` is the number of polarisations and
`<nsub>` is the number of time subintegrations.

For example `J0855-3331_2020-07-11-11:08:23_zap_chopped.16ch1p1t.ar` is a decimated archive with the edge channels removed, 16 frequency channels, 1 polarisation and 1 subintegration.


### Timing

The timing directory contains subdirectories for each project (e.g PTA) (`${params.outdir}/<pulsar>/<utc>/<beam>/decimated/<project>/`).
In each of these project timing subdirectories there is a ephemeris file (`<pulsar>.par`) and a template file (`<pulsar>.std`) which are used to make a time of arrival (ToA) `.tim` files.
The ToA files have the following format:

```
<pulsar>_<utc>_zap<_chopped>.<nchan>ch<npol>p<nsub>t.tim
```

where, as above, `<_chopped>` is included in the file name if the edge frequency channels are removed,
`<nchan>` is the number of frequency channels,
`<npol>` is the number of polarisations and
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
