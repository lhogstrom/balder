// main.nf

// Define parameters
params.processed_data_dir = '../../data/'
params.results_dir = '../../output/actionability_db_curration_20240712'
params.dbDir = '../../data/processed/balderResultsDb'
params.scripts_dir = '../code'
params.scripts1 = file('../code/clinical_annotation_db_curration.r')

log.info """\
    B A L D E R - N F   P I P E L I N E
    ===================================
    data dir     : ${params.processed_data_dir}
    database dir : ${params.dbDir}
    outdir       : ${params.results_dir}
    reports dir  : ${params.reports_dir}
    """
    .stripIndent()


// Define the preprocessing
process preprocess_data1 {
    input:
    //path raw_file from "${params.processed_data_dir}/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
    path raw_file 

    output:
    //path processed_file into processed_data1
    path "${params.dbDir}/balder-harmonized-biomarker-data-v20240712.sqlite"

    script:
    //Rscript ${params.scripts_dir}/clinical_annotation_db_curration.r $raw_file ${params.dbDir}/balder-harmonized-biomarker-data-v20240712.sqlite
    """
    Rscript $params.scripts1 $raw_file ${params.dbDir}/balder-harmonized-biomarker-data-v20240712.sqlite
    """
}

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
    preprocess_data1(file("${params.processed_data_dir}/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"))
    //preprocess_data2()
    //analysis_script1()
    //analysis_script2()
}
