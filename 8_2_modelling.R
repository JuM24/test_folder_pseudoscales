source('helper_functions.R')

# select the period for which the cumulative drug burden is to be calculated

years_used <- c(2015, 2015)
#years_used <- c(2004, 2005, 2006)

# vector with names of all four sampling appraoches
versions_all <- c('across_all', 'across_achb', 'within_all', 'within_achb')

### combine_scales ###
### if the data with the scales was exported in several chunks, combine them into a single file
for (current_version in versions_all){
  combine_scales(version = current_version, j_max = 20,  years = years_used)
}



### prepare_scales ###
### prepare the data frame for analysis (remove cases before time 0, add covariates, etc.)
for (current_version in versions_all){
  death <- prepare_scales(version = current_version, 'death', years_used)
  rm(death); gc()
  dementia <- prepare_scales(version = current_version, 'dementia', years_used)
  rm(dementia); gc()
  delirium <- prepare_scales(version = current_version, 'delirium', years_used)
  rm(delirium); gc()
}



### outcome_effect_parallel ###
### Running logistic regression models for each scale to predict the outcome using the burden score

### WARNING: this is very memory-hungry when run on several cores; 
### the `core_number` used below was the maximum number that could be stably run
### on a machine with 96 Gb RAM
death_all <- outcome_effect_parallel(version = 'across_all', 
                                     outcome_name = 'death',
                                     control = 'basic',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(), 
                                     output_file_name = 'across_all_death_unadjusted.Rds', 
                                     core_number = 10)
rm(death_all); gc()
dementia_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'dementia',
                                        control = 'basic',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_dementia_unadjusted.Rds', 
                                        core_number = 10)
rm(dementia_all); gc()
delirium_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'delirium',
                                        control = 'basic',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_delirium_unadjusted.Rds', 
                                        core_number = 10)
rm(delirium_all); gc()

## Repeat the above with adjustment
death_all <- outcome_effect_parallel(version = 'across_all', 
                                     outcome_name = 'death',
                                     control = 'full',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(), 
                                     output_file_name = 'across_all_death_adjusted.Rds', 
                                     core_number = 10)
rm(death_all); gc()
dementia_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'dementia',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_dementia_adjusted.Rds', 
                                        core_number = 10)
rm(dementia_all); gc()
delirium_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'delirium',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_delirium_adjusted.Rds', 
                                        core_number = 10)
rm(delirium_all); gc()




## modelling the effects of scales sampled "across achb"; i.e., by sampling from 
## anticholinergic drugs prescribed in the period of interest

death_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                     outcome_name = 'death',
                                     control = 'basic',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(), 
                                     output_file_name = 'across_achb_death_unadjusted.Rds', 
                                     core_number = 10)
rm(death_achb); gc()
dementia_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'dementia',
                                        control = 'basic',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_dementia_unadjusted.Rds', 
                                        core_number = 10)
rm(dementia_achb); gc()
delirium_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'delirium',
                                        control = 'basic',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_delirium_unadjusted.Rds', 
                                        core_number = 10)
rm(delirium_achb); gc()

## Repeat the above with adjustment
death_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                     outcome_name = 'death',
                                     control = 'full',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(), 
                                     output_file_name = 'across_achb_death_adjusted.Rds', 
                                     core_number = 10)
rm(death_achb); gc()
dementia_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'dementia',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_dementia_adjusted.Rds', 
                                        core_number = 10)
rm(dementia_achb); gc()
delirium_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'delirium',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_delirium_adjusted.Rds', 
                                        core_number = 10)
rm(delirium_achb); gc()






## modelling the effects of scales sampled "within all"; i.e., by sampling from 
## all drugs prescribed in the period of interest separately for each scale
death_all <- outcome_effect_parallel(version = 'within_all',
                                     outcome_name = 'death',
                                     control = 'full',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(),
                                     output_file_name = 'within_all_death_adjusted.Rds', 
                                     core_number = 5)
rm(death_all); gc()
dementia_all <- outcome_effect_parallel(version = 'within_all', 
                                        outcome_name = 'dementia',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(),
                                        output_file_name = 'within_all_dementia_adjusted.Rds', 
                                        core_number = 5)
rm(dementia_all); gc()
delirium_all <- outcome_effect_parallel(version = 'within_all', 
                                        outcome_name = 'delirium',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(),
                                        output_file_name = 'within_all_delirium_adjusted.Rds', 
                                        core_number = 4)
rm(delirium_all); gc()




## modelling the effects of scales sampled "within achb"; i.e., by sampling from 
## anticholinergic drugs prescribed in the period of interest separately for each scale
death_achb <- outcome_effect_parachbel(version = 'within_achb',
                                     outcome_name = 'death',
                                     control = 'full',
                                     smote = TRUE,
                                     model_type = 'logistic',
                                     file_path = getwd(),
                                     output_file_name = 'within_achb_death_adjusted.Rds', 
                                     core_number = 5)
rm(death_achb); gc()
dementia_achb <- outcome_effect_parachbel(version = 'within_achb', 
                                        outcome_name = 'dementia',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(),
                                        output_file_name = 'within_achb_dementia_adjusted.Rds', 
                                        core_number = 5)
rm(dementia_achb); gc()
delirium_achb <- outcome_effect_parachbel(version = 'within_achb', 
                                        outcome_name = 'delirium',
                                        control = 'full',
                                        smote = TRUE,
                                        model_type = 'logistic',
                                        file_path = getwd(),
                                        output_file_name = 'within_achb_delirium_adjusted.Rds', 
                                        core_number = 4)
rm(delirium_achb); gc()




### Running Cox proportional hazards regression for each scale to predict the outcome using the burden score

# across all
death_all <- outcome_effect_parallel(version = 'across_all', 
                                     outcome_name = 'death',
                                     control = 'basic',
                                     smote = FALSE,
                                     model_type = 'survival',
                                     competing_death = TRUE,
                                     file_path = getwd(), 
                                     output_file_name = 'across_all_death_unadjusted_Cox.Rds', 
                                     core_number = 10)
rm(death_all); gc()
dementia_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'dementia',
                                        control = 'basic',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = TRUE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_dementia_unadjusted_Cox.Rds', 
                                        core_number = 10)
rm(dementia_all); gc()
delirium_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'delirium',
                                        control = 'basic',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = TRUE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_delirium_unadjusted_Cox.Rds', 
                                        core_number = 10)
rm(delirium_all); gc()


# across achb
death_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                     outcome_name = 'death',
                                     control = 'basic',
                                     smote = FALSE,
                                     model_type = 'survival',
                                     competing_death = TRUE,
                                     file_path = getwd(), 
                                     output_file_name = 'across_achb_death_unadjusted_Cox.Rds', 
                                     core_number = 10)
rm(death_achb); gc()
dementia_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'dementia',
                                        control = 'basic',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = TRUE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_dementia_unadjusted_Cox.Rds', 
                                        core_number = 10)
rm(dementia_achb); gc()
delirium_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                        outcome_name = 'delirium',
                                        control = 'basic',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = TRUE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_achb_delirium_unadjusted_Cox.Rds', 
                                        core_number = 10)
rm(delirium_achb); gc()


## repeat above with adjusment
# across all
death_all <- outcome_effect_parallel(version = 'across_all', 
                                     outcome_name = 'death',
                                     control = 'full',
                                     smote = FALSE,
                                     model_type = 'survival',
                                     competing_death = FALSE,
                                     file_path = getwd(), 
                                     output_file_name = 'across_all_death_adjusted_Cox.Rds', 
                                     core_number = 10)
rm(death_all); gc()
dementia_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'dementia',
                                        control = 'full',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = FALSE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_dementia_adjusted_Cox.Rds', 
                                        core_number = 10)
rm(dementia_all); gc()
delirium_all <- outcome_effect_parallel(version = 'across_all', 
                                        outcome_name = 'delirium',
                                        control = 'full',
                                        smote = FALSE,
                                        model_type = 'survival',
                                        competing_death = FALSE,
                                        file_path = getwd(), 
                                        output_file_name = 'across_all_delirium_adjusted_Cox.Rds', 
                                        core_number = 10)
rm(delirium_all); gc()


# across achb
death_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                       outcome_name = 'death',
                                       control = 'full',
                                       smote = FALSE,
                                       model_type = 'survival',
                                       competing_death = FALSE,
                                       file_path = getwd(), 
                                       output_file_name = 'across_achb_death_adjusted_Cox.Rds', 
                                       core_number = 10)
rm(death_achb); gc()
dementia_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                          outcome_name = 'dementia',
                                          control = 'full',
                                          smote = FALSE,
                                          model_type = 'survival',
                                          competing_death = FALSE,
                                          file_path = getwd(), 
                                          output_file_name = 'across_achb_dementia_adjusted_Cox.Rds', 
                                          core_number = 10)
rm(dementia_achb); gc()
delirium_achb <- outcome_effect_parachbel(version = 'across_achb', 
                                          outcome_name = 'delirium',
                                          control = 'full',
                                          smote = FALSE,
                                          model_type = 'survival',
                                          competing_death = FALSE,
                                          file_path = getwd(), 
                                          output_file_name = 'across_achb_delirium_adjusted_Cox.Rds', 
                                          core_number = 10)
rm(delirium_achb); gc()