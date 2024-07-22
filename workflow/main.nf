// main.nf

// Define parameters
params.processed_data_dir = '../../data/'
params.base_dir = './'
params.results_dir = '../../output/actionability_db_curration_20240712'
params.oncokbToken = file('../../data/ONCOKB_token.txt')
params.dbDir = '../../data/processed/balderResultsDb'
params.dbName = "balder-harmonized-biomarker-data-v20240712.sqlite"
params.dbName2 = "balder-compiled-raw-data-v20240712.sqlite"
params.tmpVarTbl = "compiled_mutations_column_subset_all_studies.txt"
params.oncokbOut = "compiled_mutations_column_subset_all_studies.oncokb.txt"
params.scripts_dir = '../code'
params.scripts1 = file('../code/clinical_annotation_db_curration.r')
params.scripts2 = file('../code/cancer_genomics_studies_db_curration.r')
params.scripts3 = file('../code/setup_and_structure_oncokb_run.r')
params.scripts4 = file('../code/run_oncokb_annotator.sh')
params.scripts5 = file('../code/load_oncokb_run.r')

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

    input:
    path base_dir

    output:
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

process oncokb_file_setup {
    publishDir params.results_dir, mode:'copy'

    input:
    path base_dir
    path "$params.dbName"

    output:
    path "$params.tmpVarTbl"

    script:
    """
    Rscript $params.scripts3 $base_dir $params.dbName $params.tmpVarTbl
    """
}

process run_oncokb {
    publishDir params.results_dir, mode:'copy'

    input:
    path "$params.tmpVarTbl"
    path "$params.oncokbToken"

    output:
    path "$params.oncokbOut"

    script:
    """
    bash $params.scripts4 $params.tmpVarTbl $params.oncokbToken $params.oncokbOut
    """
}

process load_oncokb_output {
    publishDir params.results_dir, mode:'copy'
    
    input:
    path base_dir
    path "$params.dbName"

    output:
    path "$params.dbName"

    script:
    """
    Rscript $params.scripts5 $base_dir $params.dbName
    """
}

workflow {
    sqldb = preprocess_data1(file("${params.base_dir}"))
    sqldb2 = preprocess_data2(file("${params.base_dir}"),sqldb)
    tmpVarTbl = oncokb_file_setup(file("${params.base_dir}"),sqldb2)
    //oncokbOut = run_oncokb(tmpVarTbl,file("${params.oncokbToken}"))
    sqldb3 = load_oncokb_output(file("${params.base_dir}"),sqldb2)
}

