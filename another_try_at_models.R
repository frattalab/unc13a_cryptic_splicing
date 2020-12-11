library(tidymodels)
library(usemodels)

set.seed(42)
clean_data_table = data.table::fread(file.path(here::here(),"data","nygc_junction_information.csv"))

model_df <- clean_data_table %>% 
    filter(UNC13A_annotated_leaf > 10) %>% 
    filter(grepl("S-TDP|D-TDP",disease_group2)) %>% 
    filter(disease_tissue == T) %>% 
    mutate(stmn2_psi_by_number_g_alleles = (number_g_alleles + 1) * stmn_2_cryptic_psi_leaf) %>% 
    dplyr::select(age,
                  unc13a_cryptic_leaf_psi,
                  stmn2_psi_by_number_g_alleles,
                  stmn_2_cryptic_psi_leaf,
                  number_g_alleles,
                  disease_group2,
                  tissue_clean,
                  reported_mutations,
                  sex) %>% 
    mutate(disease_group2 = snakecase::to_snake_case(disease_group2)) %>% 
    dplyr::rename(disease_group = disease_group2)



model_split <- initial_split(model_df, strata = unc13a_cryptic_leaf_psi)
model_train <- training(model_split)
model_test <- testing(model_split)

set.seed(234)
model_folds <- bootstraps(model_train, strata = unc13a_cryptic_leaf_psi)

# now a random forest -----------------------------------------------------

ranger_recipe <- 
    recipe(formula = unc13a_cryptic_leaf_psi ~ ., data = model_df) %>% 
    step_string2factor(one_of("disease_group", "tissue_clean", "reported_mutations", 
                              "sex")) 

ranger_spec <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
    set_mode("regression") %>% 
    set_engine("ranger") 

ranger_workflow <- 
    workflow() %>% 
    add_recipe(ranger_recipe) %>% 
    add_model(ranger_spec) 

set.seed(10356)
doParallel::registerDoParallel()
ranger_tune <-
    tune_grid(ranger_workflow, resamples = model_folds,
              grid = 10)
show_best(ranger_tune, metric = "rmse")
show_best(ranger_tune, metric = "rsq")

autoplot(ranger_tune)

final_ranger <- ranger_workflow %>% 
    finalize_workflow(select_best(ranger_tune, "rmse"))

ranger_last <- last_fit(final_ranger, model_split)
collect_metrics(ranger_last)

collect_predictions(ranger_last) %>% 
    ggplot(aes(x = unc13a_cryptic_leaf_psi, y = .pred)) + 
    geom_abline() + 
    geom_point(pch = 21) + 
    coord_fixed()

test_ex <- bake(prep(ranger_recipe), new_data = model_df)
x_data <- test_ex %>% dplyr::select(-unc13a_cryptic_leaf_psi) %>% as.data.frame()

ranger_predictor <- iml::Predictor$new(fit(final_ranger, model_df) %>% pull_workflow_fit(), data = x_data, y = test_ex$unc13a_cryptic_leaf_psi)
ranger_feature_imp <- iml::FeatureImp$new(ranger_predictor, loss = "rmse")
ranger_feature_imp %>% plot() + xlim(1,1.3) + xlab("Feature Importance (RMSE)")
ranger_ale <- FeatureEffects$new(ranger_predictor)
ranger_ale %>% plot() +
    scale_x_discrete(guide = guide_axis(n.dodge=2))
