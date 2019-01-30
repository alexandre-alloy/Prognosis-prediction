#!/bin/bash
RSCRIPT=/nfs/apps/R/3.5.1/bin/Rscript
FUNCTIONS_RDATA=/ifs/scratch/c2b2/ac_lab/apa2118/rscripts/classification_by_MRs/Functions.RData #RData file with all functions
LIST_PATIENT_GROUPS=/ifs/scratch/c2b2/ac_lab/apa2118/rscripts/classification_by_MRs/AML/list_patients.rds #rda file with a list of lists of patients. use patient categories as names
PATIENT_GROUP_1="multiple_intersect(High_risk, diagnostic, AAML0531)" #testing group, group with a phenotype, e.g. "dead", "relapse"
PATIENT_GROUP_2="multiple_intersect(Low_risk, diagnostic, AAML0531)" #name of a list corresponding to the first group of patients (control group, e.g "alive", "censored")
EXP_NAME=AML_HR_LR_diag_AAML0531
WORKING_DIR=/ifs/scratch/c2b2/ac_lab/apa2118/rscripts/classification_by_MRs/AML/${EXP_NAME}
PATH_TO_GENE_EXP_MAT=/ifs/scratch/c2b2/ac_lab/apa2118/rscripts/classification_by_MRs/AML/AML_cpm_entrez.rds #in txt or rds format
PATH_TO_INTERACTOME=/ifs/scratch/c2b2/ac_lab/apa2118/rscripts/classification_by_MRs/AML/interactome_TF_coTF_AML.rds #in rds format
NUMBER_OF_FOLDS=5
NUMBER_OF_CLUSTERS=1
CLUSTERING_ALGORITHM=Mclust #available options: consensus clustering, Mclust, none (if NUMBER_OF_CLUSTERS == 1: "none" is assigned automatically)
CLUSTER_GEM=TRUE #if TRUE, clustering is performed on the GEM, if FALSE clustering is performed on the VIPER matrix
CLASSIFICATION_ALGORITHM=lda #available options: random_forest, logistic_regression, lda, svm
PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION=0.25
CV_TYPE=monte_carlo #available options: kfold , LOOCV, monte_carlo
VALIDATION_FOLD_SIZE=0.25 #Only used with monte carlo - a proportion of the sample set for validation e.g. 0.3
COMPUTE_NULL=TRUE
RANDOM_NEGATIVE_CONTROL=TRUE #if false, "non-significant" features are used as predictors for the negative control
EQUILIBRATE_CLASSES=FALSE
DOWNSAMPLE_MAJORITY_CLASS=FALSE
OVERSAMPLE_MINORITY_CLASS=FALSE
SMOTE=FALSE
TOMEK=FALSE
RUN_VST=FALSE
mRNA_CONTROL=FALSE
TOP_AND_BOTTOM_FEATURES=FALSE  #if FALSE, will take top features only
EQUILIBRATE_TOP_AND_BOTTOM_FEATURES=FALSE #if TRUE, will take the same number of top and bottom features
NB_PIPELINE_ITERATIONS=$(bc<<<"1/${PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION}")
#remainder=$(bc<<<"1%${PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION}")
# if [  bc <<< "(1%${PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION})!=0" ]; then
#   echo "ERROR: variable PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION must divide 1 with no remainder, e.g. 0.1, 0.2, 0.25, etc"
#   exit
# fi


########################################################
mkdir -p ${WORKING_DIR}

#STEP 1: assign patients to folds
echo "load('${FUNCTIONS_RDATA}')
if ('txt' %in% unlist(strsplit('${PATH_TO_GENE_EXP_MAT}', split = '\\\.'))) GEM = read.delim(file='${PATH_TO_GENE_EXP_MAT}', header=T, row.names=1, sep='\t', as.is = T)
if ('rds' %in% unlist(strsplit('${PATH_TO_GENE_EXP_MAT}', split = '\\\.'))) GEM = readRDS('${PATH_TO_GENE_EXP_MAT}')

list_patients = readRDS('${LIST_PATIENT_GROUPS}')
for (i in names(list_patients)) assign(x = as.character(i), value = list_patients[[i]])
patient_group_1 = ${PATIENT_GROUP_1}
patient_group_2 = ${PATIENT_GROUP_2}

Interactome <- readRDS('${PATH_TO_INTERACTOME}')

Initial_processing_folds(patient_group_1 = patient_group_1,
                         patient_group_2 = patient_group_2,
                         CV_type = '${CV_TYPE}',
                         nb_folds = ${NUMBER_OF_FOLDS},
                         validation_fold_size = ${VALIDATION_FOLD_SIZE},
                         holdout_data = ${PCT_HELDOUT_DATA_PER_PIPELINE_ITERATION},
                         equilibrate_classes = ${EQUILIBRATE_CLASSES},
                         downsample_majority_class = ${DOWNSAMPLE_MAJORITY_CLASS},
                         oversample_minority_class = ${OVERSAMPLE_MINORITY_CLASS},
                         SMOTE = ${SMOTE},
                         TOMEK = ${TOMEK},
                         GEM = GEM,
                         interactome = Interactome,
                         save_path = './intermediate_files',
                         savefile = TRUE)
      " > ${WORKING_DIR}/Initial_processing_folds.r

 chmod 755 ${WORKING_DIR}/Initial_processing_folds.r

 qsub  -N Initial_processing_folds_${EXP_NAME} \
       -l mem=3G,time=1:: \
       -wd ${WORKING_DIR} \
       -b y ${RSCRIPT} \
       ${WORKING_DIR}/Initial_processing_folds.r

#STEP 2: normalization and clustering
for ((pipeline_iter=1;pipeline_iter<=${NB_PIPELINE_ITERATIONS};pipeline_iter++))
do
  for ((folds_iter=1;folds_iter<=${NUMBER_OF_FOLDS};folds_iter++))
  do
    echo "load('${FUNCTIONS_RDATA}')
    Initial_normalization_clust(Initial_processing_folds_object,
                                        GEM = NULL,
                                        interactome = NULL,
                                        compute_folds = ${folds_iter},
                                        nb_clusters = ${NUMBER_OF_CLUSTERS},
                                        cluster_GEM = ${CLUSTER_GEM},
                                        run_vst = ${RUN_VST},
                                        clustering_algo = '${CLUSTERING_ALGORITHM}',
                                        clustering_subset = NULL,
                                        mRNA_control = ${mRNA_CONTROL},
                                        pipeline_iter = ${pipeline_iter},
                                        save_path = './intermediate_files',
                                        savefile = TRUE,
                                        loadfile = TRUE)
" > ${WORKING_DIR}/initial_normalization_clust_p${pipeline_iter}_k${folds_iter}.r

    chmod 755 ${WORKING_DIR}/initial_normalization_clust_p${pipeline_iter}_k${folds_iter}.r

    qsub -N initial_normalization_clust_p${pipeline_iter}_k${folds_iter}_${EXP_NAME} \
         -hold_jid Initial_processing_folds_${EXP_NAME} \
         -l mem=6G,time=4:: \
         -wd ${WORKING_DIR} \
         -b y ${RSCRIPT} \
         ${WORKING_DIR}/initial_normalization_clust_p${pipeline_iter}_k${folds_iter}.r
  done
done


#Step 3: validation
for ((pipeline_iter=1;pipeline_iter<=${NB_PIPELINE_ITERATIONS};pipeline_iter++))
do
  for ((folds_iter=1;folds_iter<=${NUMBER_OF_FOLDS};folds_iter++))
  do
    for ((clust_iter=1;clust_iter<=${NUMBER_OF_CLUSTERS};clust_iter++))
    do
      echo "load('${FUNCTIONS_RDATA}')
            validation(Risk_group_classification_object = NULL,
                        k = ${folds_iter},
                        c = ${clust_iter},
                        pipeline_iter = ${pipeline_iter},
                        classification_algorithm = '${CLASSIFICATION_ALGORITHM}',
                        nb_iterations = 100,
                        MR_range = 1:50,
                        interactome = NULL,
                        GEM = NULL,
                        top_and_bottom = ${TOP_AND_BOTTOM_FEATURES},
                        equilibrate_top_bottom = ${EQUILIBRATE_TOP_AND_BOTTOM_FEATURES},
                        mRNA_control = ${mRNA_CONTROL},
                        random_negative_control = ${RANDOM_NEGATIVE_CONTROL},
                        compute_null = ${COMPUTE_NULL},
                        equilibrate_classes = FALSE, #different from folds- equilibrate classes. Here will drop off extra patients from the majority class
                        save_path = './intermediate_files',
                        savefile = TRUE,
                        loadfile = TRUE)" > ${WORKING_DIR}/validation_p${pipeline_iter}_k${folds_iter}_c${clust_iter}_${EXP_NAME}_${CLASSIFICATION_ALGORITHM}.r

            chmod 755 ${WORKING_DIR}/validation_p${pipeline_iter}_k${folds_iter}_c${clust_iter}_${EXP_NAME}_${CLASSIFICATION_ALGORITHM}.r

            qsub -N validation_p${pipeline_iter}_k${folds_iter}_c${clust_iter}_${EXP_NAME}_${CLASSIFICATION_ALGORITHM} \
            -hold_jid initial_normalization_clust_p${pipeline_iter}_k${folds_iter}_${EXP_NAME} \
            -l mem=12G,time=4:: \
            -wd ${WORKING_DIR} \
            -b y ${RSCRIPT} ${WORKING_DIR}/validation_p${pipeline_iter}_k${folds_iter}_c${clust_iter}_${EXP_NAME}_${CLASSIFICATION_ALGORITHM}.r
    done
  done
done

#Step 4: consolidate all results and create final list of features for testing
echo "load('${FUNCTIONS_RDATA}')
consolidate <- merge_norm_clust_validation_files()
MRs_per_fold(Risk_group_classification_object = consolidate,
                         k = 1:"${NUMBER_OF_FOLDS}",
                         c = 1:"${NUMBER_OF_CLUSTERS}",
                         pipeline_iter = 1:"${NB_PIPELINE_ITERATIONS}",
                         p_value_threshold = 0.05,
                         use_null = ${COMPUTE_NULL},
                         savefile = TRUE,
                         loadfile = TRUE)
" > ${WORKING_DIR}/feature_consolidation.r

chmod 755 ${WORKING_DIR}/feature_consolidation.r

qsub -N feature_consolidation_${EXP_NAME}_${CLASSIFICATION_ALGORITHM} \
-hold_jid validation_*${EXP_NAME}_${CLASSIFICATION_ALGORITHM} \
-l mem=6G,time=1:: \
-wd ${WORKING_DIR} \
-b y ${RSCRIPT} ${WORKING_DIR}/feature_consolidation.r

#Step 5: final validation in test set
for ((pipeline_iter=1;pipeline_iter<=${NB_PIPELINE_ITERATIONS};pipeline_iter++))
do
  for ((clust_iter=1;clust_iter<=${NUMBER_OF_CLUSTERS};clust_iter++))
  do
    echo "load('${FUNCTIONS_RDATA}')
    Internal_validation(Risk_group_classification_object = NULL,
                                pipeline_iter = ${pipeline_iter},
                                nb_folds = ${NUMBER_OF_FOLDS},
                                c = ${NUMBER_OF_CLUSTERS},
                                classification_algorithm = '${CLASSIFICATION_ALGORITHM}',
                                nb_iterations = 100,
                                interactome,
                                mRNA_control = ${mRNA_CONTROL},
                                compute_null = ${COMPUTE_NULL},
                                savefile = TRUE,
                                loadfile = TRUE)" > ${WORKING_DIR}/internal_validation_p${pipeline_iter}_c${clust_iter}.r

     chmod 755 ${WORKING_DIR}/internal_validation_p${pipeline_iter}_c${clust_iter}.r

     qsub -N internal_validation_p${pipeline_iter}_c${clust_iter}_${EXP_NAME} \
     -hold_jid feature_consolidation_${EXP_NAME}_${CLASSIFICATION_ALGORITHM} \
     -l mem=6G,time=4:: \
     -wd ${WORKING_DIR} \
     -b y ${RSCRIPT} ${WORKING_DIR}/internal_validation_p${pipeline_iter}_c${clust_iter}.r

  done
done

#Step 6: final consolidation
echo "load('${FUNCTIONS_RDATA}')
Final_consolidation(results_file_name = '${EXP_NAME}.rda')" > ${WORKING_DIR}/final_result.r

chmod 755 ${WORKING_DIR}/final_result.r

qsub -N final_result_${EXP_NAME} \
-hold_jid internal_validation*${EXP_NAME} \
-l mem=3G,time=1:: \
-wd ${WORKING_DIR} \
-b y ${RSCRIPT} ${WORKING_DIR}/final_result.r
