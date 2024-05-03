#!/usr/bin/env bash
ONCKBPATH=/Users/larsonhogstrom/Documents/code/oncokb-annotator

#subset
#IMAF="../../output/clinical_annotation_matching_20240412/scratch_compiled_mutations_sample_info_all_studies.txt"
#OMAF="../../output/clinical_annotation_matching_20240412/scratch_compiled_mutations_sample_info_all_studies.oncokb.txt"

# subset of columns, all HGVSp variants
#IMAF="../../output/clinical_annotation_matching_20240412/compiled_mutations_column_subset_all_studies.txt"
#OMAF="../../output/clinical_annotation_matching_20240412/compiled_mutations_column_subset_all_studies.oncokb.txt"

# scratch subset of lines for debugging
IMAF="../../output/clinical_annotation_matching_20240412/compiled_mutations_column_subset_all_studies_columns_979400_to_2174415.txt"
OMAF="../../output/clinical_annotation_matching_20240412/compiled_mutations_column_subset_all_studies_columns_979400_to_2174415.oncokb.txt"


#IC="data/example_clinical.txt"
#OC="data/example_clinical.oncokb.txt"

TOKEN="06693fc8-bbb2-4b75-ae0a-728fa96be6e9" #OncoKB API Token
README="data/example_README.txt"

#python $ONCKBPATH/MafAnnotator.py -i "$IMAF" -o "$OMAF" -c "$IC" -b "$TOKEN"
python $ONCKBPATH/MafAnnotator.py -i "$IMAF" -o "$OMAF" -b "$TOKEN"

#python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSPSHORT" -c "$IC" -b "$TOKEN" -q hgvsp_short
#python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSP" -c "$IC" -b "$TOKEN" -q hgvsp
#python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSG" -c "$IC" -b "$TOKEN" -q hgvsg
#python MafAnnotator.py -i "$IMAF" -o "$OMAFGC" -c "$IC" -b "$TOKEN" -q genomic_change
#
#python MafAnnotator.py -i "$IMAF38" -o "$OMAF38" -c "$IC" -b "$TOKEN"
#
#python MafAnnotator.py -i "$IATYPICALALT" -o "$OATYPICALALT" -c "$IC" -b "$TOKEN"
#
#python FusionAnnotator.py -i "$IF" -o "$OF" -c "$IC" -b "$TOKEN"
#python StructuralVariantAnnotator.py -i "$ISV" -o "$OSV" -c "$IC" -b "$TOKEN"
#python CnaAnnotator.py -i "$ICNA" -o "$OCNA" -c "$IC" -b "$TOKEN"
#python CnaAnnotator.py -i "$IICNA" -o "$OICNA" -c "$IC" -b "$TOKEN" -f "individual" -z
#python ClinicalDataAnnotator.py -i "$IC" -o "$OC" -a "$OMAF,$OATYPICALALT,$OCNA,$OF,$OSV"
#
#python GenerateReadMe.py -o "$README"
