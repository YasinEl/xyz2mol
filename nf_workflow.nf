#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.download_link = '' // No default download link
params.csv_root_dir = '' // No default csv root dir
params.xyz_root_dir = 'TMPQCXMS'

TOOL_FOLDER = "$baseDir"

process DownloadData {
    when:
    params.download_link != ''

    output:
    path("tar_files/*"), emit: tarFile

    script:
    """
    mkdir -p tar_files
    wget -P tar_files ${params.download_link}
    """
}

process ExtractTar {
    input:
    path tar_file 

    output:
    path("${params.xyz_root_dir}"), emit: extractedFiles

    script:
    """
    tar -xf ${tar_file} -C .
    """
}


process RenameFiles {
    input:
    path inputDir
    val toolFolder

    output:
    path("${inputDir}/*"), emit: renamedFiles

    script:
    """
    python3 $toolFolder/rename_files.py -directory "$inputDir"
    """
}


process ConvertXYZtoCSV {
    input:
    path xyz_file 
    val toolFolder

    output:
    path "*.csv", emit: csvFile

    script:
    """
    csv_file="\${xyz_file/unique_named_xyzs/params.csv_root_dir}"
    csv_file="\${csv_file%.xyz}.csv"
    csv_dir=\$(dirname "\$csv_file")
    mkdir -p "\$csv_dir"
    python3 $toolFolder/main.py -xyz "\$xyz_file" -csv "\$csv_file"
    """
}

workflow {
    tarFiles = DownloadData()
    extractedFiles = ExtractTar(tarFiles)
    renamedFiles = RenameFiles(extractedFiles, TOOL_FOLDER)
    ConvertXYZtoCSV(renamedFiles, TOOL_FOLDER)
}

