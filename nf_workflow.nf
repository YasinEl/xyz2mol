#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.download_link = '' // No default download link
params.csv_root_dir = './csvs' // No default csv root dir
params.out_root_dir = './out' // No default csv root dir
params.xyz_root_dir = 'TMPQCXMS'
params.download_links = ""
download_links_list = params.download_links.split(',')
download_links_channel = Channel.from(download_links_list)
download_links_channel.view()

TOOL_FOLDER = "$baseDir"

process DownloadData {
    input:
    each link

    output:
    path("tar_files/*"), emit: tarFile

    script:
    """
    mkdir -p tar_files
    wget -P tar_files $link
    """
}

process ExtractTar {
    input:
    each tar_files

    output:
    path "TMPQCXMS_${tar_files.baseName}/*", emit: extractedFiles
    // path "TMPQCXMS_*/*", emit: extractedFiles

    """
    tar -xf ${tar_files} -C .
    mv ./TMPQCXMS ./TMPQCXMS_${tar_files.baseName}
    """
}



process RenameFiles {
    cache 'lenient'

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
    cache 'lenient'

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
    cache 'lenient'
    
    conda "$TOOL_FOLDER/requirements.yml"

    publishDir "./csvs", mode: 'copy'
    //errorStrategy 'ignore'

    input:
    each xyz_file 
    val toolFolder

    output:
    path "*.csv", emit: csvFile, optional: true

    script:
    """
    csv_file="${xyz_file.baseName}.csv"


    # Check if the file has a .xyz extension
    if [[ "$xyz_file" == *.xyz  && "$xyz_file" != *CID* && "$xyz_file" != *MDtrj* ]]; then
        python3 $toolFolder/main.py -xyz "$xyz_file" -csv "\$csv_file"
    else
        echo "Not an .xyz file. Skipping..."
    fi
    """
}

process ParseOutFiles {
    cache 'lenient'

    conda "$TOOL_FOLDER/requirements.yml"

    // publishDir "./out", mode: 'symlink'

    input:
    // path inputDir  // Directly specify the directory
    val toolFolder
    each outfile

    output:
    path "*.json", emit: jsonFile, optional: true

    script:
    """
    if [[ "$outfile" == *.out ]]; then
        python3 $toolFolder/parse_out_file.py --filepath "$outfile"
    else
        echo "Not an .out file. Skipping..."
    fi
    """
}

process SummarizeTrajectories {
    cache 'lenient'

    conda "$TOOL_FOLDER/requirements.yml"

    publishDir "./summary_csvs", mode: 'symlink'

    input:
    // path inputDir  // Directly specify the directory
    val toolFolder
    val csvfiles

    output:
    path "*.csv", emit: csvFile, optional: true

    script:
    """
    mkdir -p summary_csvs
    python3 $toolFolder/summarize_trajectories.py --input "${workflow.launchDir}/${params.csv_root_dir}"
    """
}

process SummarizeSinglets {
    cache 'lenient'

    conda "$TOOL_FOLDER/requirements.yml"

    publishDir "./singlet_csvs", mode: 'copy' 

    input:
    // path inputDir  // Directly specify the directory
    val toolFolder
    val csvfiles
    each outFiles

    output:
    path "*.csv", emit: csvFile, optional: true

    script:
    """
    python3 $toolFolder/summarize_singlets.py --input_json "$outFiles" --input_csvs "${workflow.launchDir}/${params.csv_root_dir}"
    """
}

process CreateNetworkBase {
    cache 'lenient'

    conda "$TOOL_FOLDER/requirements.yml"

    publishDir "./summary_csvs", mode: 'copy' 

    input:
    // path inputDir  // Directly specify the directory
    val toolFolder
    val csvfiles

    output:
    path "*.csv", emit: csvFile, optional: true

    script:
    """
    python3 $toolFolder/summarize_network.py --input_path "${workflow.launchDir}/singlet_csvs"
    """
}

 workflow {
     download_links_channel
     tarFiles = DownloadData(download_links_channel)
     ExtractTar(tarFiles)
     RenameFiles(ExtractTar.out.extractedFiles.collect().map{ it.sort{ a, b -> a.toString() <=> b.toString() } }.flatten(), TOOL_FOLDER)
     renamed_files = AddDirectoryToNames(RenameFiles.out.renamedFiles.collect().map{ it.sort{ a, b -> a.toString() <=> b.toString() } }.flatten(), TOOL_FOLDER)
     outFiles = ParseOutFiles(TOOL_FOLDER, renamed_files.collect()).jsonFile   
     csvFiles = ConvertXYZtoCSV(AddDirectoryToNames.out.renamedFiles.collect().map{ it.sort{ a, b -> a.toString() <=> b.toString() } }.flatten(), TOOL_FOLDER).csvFile
     singlets = SummarizeSinglets(TOOL_FOLDER, csvFiles.collect(), ParseOutFiles.out.jsonFile.collect()).csvFile
     network_data = CreateNetworkBase(TOOL_FOLDER, singlets.collect())
 }



