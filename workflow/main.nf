// main.nf

// Define parameters
params.processed_data_dir = '../../data/'
params.base_dir = './'
params.results_dir = '../../output/actionability_db_curration_20240712'
params.dbDir = '../../data/processed/balderResultsDb'
params.dbName = "balder-harmonized-biomarker-data-v20240712.sqlite"
params.dbName2 = "balder-compiled-raw-data-v20240712.sqlite"
params.scripts_dir = '../code'
params.scripts1 = file('../code/clinical_annotation_db_curration.r')
params.scripts2 = file('../code/cancer_genomics_studies_db_curration.r')

log.info """\
    B A L D E R - N F   P I P E L I N E
    ===================================
    data dir     : ${params.processed_data_dir}
    base dir     : ${params.base_dir}
    database dir : ${params.dbDir}
    outdir       : ${params.results_dir}
    reports dir  : ${params.reports_dir}
    """
    .stripIndent()


// Define the preprocessing
process preprocess_data1 {
    //publishDir params.dbDir, mode:'copy'

    input:
    //path raw_file from "${params.processed_data_dir}/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
    //path raw_file 
    path base_dir

    output:
    //path processed_file into processed_data1
    path "$params.dbName"

    script:
    """
    Rscript $params.scripts1 $base_dir $params.dbName
    """
}

process preprocess_data2 {
    publishDir params.dbDir, mode:'copy'

    input:
    path base_dir
    path "$params.dbName"

    output:
    path "$params.dbName"

    script:
    """
    Rscript $params.scripts2 $base_dir $params.dbName $params.dbName2
    """
}


// Define the preprocessing
//process preprocess_data1 {
//    //publishDir("$params.results_folder", mode: "copy", overwrite: false)
//
//    input:
//        val results_filepath
//        val true_proportions_path
//    output:
//        path("deconvolution_analysis_*")
//
//    script:
//    """
//    Rscript $params.analyze_results_script -r $results_filepath -t $true_proportions_path
//    """
//}


//process preprocess_data2 {
//    input:
//    path processed_file from processed_data1
//
//    output:
//    path processed_file2 into processed_data2
//
//    script:
//    """
//    Rscript ${params.scripts_dir}/preprocessing/preprocess_data2.R ${params.processed_data_dir}/processed_data1.csv ${params.processed_data_dir}/processed_data2.csv
//    """
//}
//
//// Define the analysis process
//process analysis_script1 {
//    input:
//    path processed_file2 from processed_data2
//
//    output:
//    path result1 into analysis1_result
//
//    script:
//    """
//    Rscript ${params.scripts_dir}/analysis/analysis_script1.R ${params.processed_data_dir}/processed_data2.csv ${params.results_dir}/result1.csv
//    """
//}
//
//process analysis_script2 {
//    input:
//    path processed_file2 from processed_data2
//
//    output:
//    path result2 into analysis2_result
//
//    script:
//    """
//    Rscript ${params.scripts_dir}/analysis/analysis_script2.R ${params.processed_data_dir}/processed_data2.csv ${params.results_dir}/result2.png
//    """
//}

// Workflow definition
workflow {
    //preprocess_data1(file("${params.processed_data_dir}/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"))
    sqldb = preprocess_data1(file("${params.base_dir}"))
    sqldb2 = preprocess_data2(file("${params.base_dir}"),sqldb)
    //analysis_script1()
    //analysis_script2()
}
