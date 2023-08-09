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
    path("${inputDir}"), emit: renamedFiles

    script:
    """
    python3 $toolFolder/rename_files.py -directory "$inputDir"
    """
}

process AddDirectoryToNames {
    input:
    path inputDir
    val toolFolder

    output:
    path("${inputDir}/**"), emit: renamedFiles

    script:
    """
    python3 $toolFolder/add_directory_name.py -directory "$inputDir"
    """
}

process ConvertXYZtoCSV {
    conda "$TOOL_FOLDER/requirements.yml"

    publishDir "./csvs", mode: 'copy'

    input:
    each xyz_file 
    val toolFolder

    output:
    path "*.csv", emit: csvFile, optional: true

    script:
    """
    csv_file="${xyz_file.baseName}.csv"


    # Check if the file has a .xyz extension
    if [[ "$xyz_file" == *.xyz ]]; then
        python3 $toolFolder/main.py -xyz "$xyz_file" -csv "\$csv_file"
    else
        echo "Not an .xyz file. Skipping..."
    fi

    
    
    
    """
}

workflow {
    tarFiles = DownloadData()
    extractedFiles = ExtractTar(tarFiles)
    renamedFiles = RenameFiles(extractedFiles, TOOL_FOLDER)
    directroy_Files = AddDirectoryToNames(renamedFiles, TOOL_FOLDER)
    ConvertXYZtoCSV(directroy_Files, TOOL_FOLDER)
}

