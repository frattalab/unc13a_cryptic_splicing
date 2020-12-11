library(tidymodels) # Includes the recipes package
library(tidyposterior)
library(rstanarm)
set.seed(42)
clean_data_table = fread(file.path(here::here(),"data","nygc_junction_information.csv"))

model_df <- clean_data_table %>% 
    filter(unc13a_annotated > 10) %>% 
    filter(grepl("S-TDP|D-TDP",disease_group2)) %>% 
    filter(disease_tissue == T) %>% 
    mutate(unc13a_cryptic_leaf_psi_log10 = log10(unc13a_cryptic_leaf_psi + 1)) %>% 
    mutate(number_g_alleles = number_g_alleles + 1) %>% 
    dplyr::select(sample,age,
    unc13a_cryptic_leaf_psi_log10,
    stmn_2_cryptic_psi_leaf,
    number_g_alleles,
    disease_group2,
    tissue_clean,
    reported_mutations,
    sex)

model_split <- initial_split(model_df, prob = 0.80, strata = unc13a_cryptic_leaf_psi_log10)
model_train <- training(model_split)
model_test  <-  testing(model_split)

####log transform and 
simple_rec <- 
    recipe(unc13a_cryptic_leaf_psi_log10 ~ .,
           data = model_train) %>%
    update_role(sample, new_role = "sample_id") %>% 
    step_dummy(all_predictors()) %>% 
    step_log(stmn_2_cryptic_psi_leaf, offset = 1)  %>%  
    step_interact(terms =  ~ stmn_2_cryptic_psi_leaf:number_g_alleles)

# models specs ----------------------------------------------------------

standard_lm <- linear_reg() %>% 
    set_engine("lm")

# linear workflow  ----------------------------------------------------------

lm_wf <- 
    workflow() %>% 
    add_model(standard_lm) %>% 
    add_recipe(simple_rec)

lm_fit <- fit(lm_wf, model_train)
lm_fit %>% 
    # This returns the parsnip object:
    pull_workflow_fit() %>% 
    # Now tidy the linear model object:
    tidy() %>% 
    filter(p.value < 0.1)


# random forest workflow -----------------------------------------------------------

rf_model <- 
    rand_forest(trees = 1000) %>% 
    set_engine("ranger") %>% 
    set_mode("regression")

rf_wf <- 
    workflow() %>% 
    add_formula(
        unc13a_cryptic_leaf_psi_log10 ~ .) %>% 
    add_model(rf_model) 


# xgboost workflow --------------------------------------------------------

XGBoost_model <- boost_tree() %>% 
    set_mode("regression") %>% 
    set_engine("xgboost")



# now we're doing folds ---------------------------------------------------
set.seed(42)
model_folds <- vfold_cv(model_train, v = 10,
                        repeats = 5,
                        strata = unc13a_cryptic_leaf_psi_log10)


ames_folds <- vfold_cv(model_split, v = 10,strata = unc13a_cryptic_leaf_psi_log10)

keep_pred <- control_resamples(save_pred = TRUE)

set.seed(130)
library(doMC)
registerDoMC(cores = 2)

tidy_k_folds <- vfold_cv(model_train)

XGBoost_res <- fit_resamples(XGBoost_model, simple_rec, model_folds,control = keep_pred)

rf_res <- 
    rf_wf %>% 
    fit_resamples(resamples = model_folds, control = keep_pred)

lm_res <-
    lm_wf %>%
    fit_resamples(resamples = model_folds, control = keep_pred)

collect_metrics(XGBoost_res)
collect_metrics(rf_res)
collect_metrics(lm_res)

xgb_rsq <- 
    collect_metrics(XGBoost_res, summarize = FALSE) %>% 
    filter(.metric == "rsq") %>% 
    dplyr::select(id, id2, `xgb model` = .estimate)

lm_rsq <- 
    collect_metrics(lm_res, summarize = FALSE) %>% 
    filter(.metric == "rsq") %>% 
    dplyr::select(id, id2, `linear model` = .estimate)

rf_rsq <- 
    collect_metrics(rf_res, summarize = FALSE) %>% 
    filter(.metric == "rsq") %>% 
    dplyr::select(id, id2,`random forest` = .estimate)  

rsq_estimates <- 
    inner_join(lm_rsq, rf_rsq, by = c("id","id2"))  %>% 
    inner_join(xgb_rsq,by = c("id","id2"))  %>% 
    mutate(id_rp = paste0(id,"_",id2))

rsq_estimates %>% 
    dplyr::select(-id,-id2) %>% 
    data.table::melt(id.vars = 'id_rp') %>% 
    mutate(model = reorder(variable, value)) %>% 
    ggplot(aes(x = variable, y = value, group = id_rp, col = id_rp)) + 
    geom_line(alpha = .5, lwd = 1.25) + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = expression(paste(R^2, "statistics")))

# tuning my models --------------------------------------------------------


rf_model_tunable <- 
    rand_forest(mtry = tune(),
                trees = tune(),
                min_n = tune()) %>% 
    set_engine("ranger") %>% 
    set_mode("regression")


xgb_model_tunable <- boost_tree(mtry = tune(),
                                trees = tune(),
                                tree_depth = tune(),
                                learn_rate = tune()) %>% 
    set_mode("regression") %>% 
    set_engine("xgboost")

####make the gride of parameters to try####
forest_grid <- bind_cols(tibble(mtry = rep(c(3, 4, 5,6),4)),grid_regular( trees(), min_n(), levels = 4))
xgb_grid <- bind_cols(tibble(mtry = rep(c(3, 4, 5,6,7,8),36)),grid_regular( trees(), tree_depth(),learn_rate(), levels = 6))
####tune the grid####
doParallel::registerDoParallel()
set.seed(345)
forest_rs <- tune_grid(
    rf_model_tunable,
    unc13a_cryptic_leaf_psi_log10 ~ .,
    resamples = model_folds,
    grid = forest_grid,
    metrics = metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae)
)
collect_metrics(forest_rs)
show_best(forest_rs, "rmse")
final_forest <- finalize_model(rf_model_tunable, select_best(forest_rs, "rmse"))
tuned_rf_res <- 
    workflow() %>% 
    add_formula(
        unc13a_cryptic_leaf_psi_log10 ~ .) %>% 
    add_model(final_forest)  %>% 
    fit_resamples(resamples = model_folds, control = keep_pred)

tuned_rf_rsq <- 
    collect_metrics(tuned_rf_res, summarize = FALSE) %>% 
    filter(.metric == "rsq") %>% 
    dplyr::select(id, id2,`random forest tuned` = .estimate)  

rsq_estimates <- rsq_estimates %>% 
    inner_join(tuned_rf_rsq,by = c("id","id2"))  %>% 
    mutate(id_rp = paste0(id,"_",id2))

rmse_estimates <-     collect_metrics(rf_res, summarize = FALSE) %>% 
    filter(.metric == "rmse") %>% 
    dplyr::select(id, id2,`random forest` = .estimate)  %>% 
    inner_join(
        collect_metrics(tuned_rf_res, summarize = FALSE) %>% 
            filter(.metric == "rmse") %>% 
            dplyr::select(id, id2,`tuned random forest` = .estimate),
        by = c("id","id2")
    )

# did tuning make the tree better? ----------------------------------------

compare_models_rsq <- 
    model_folds %>% 
    left_join(rsq_estimates) %>% 
    dplyr::select(-id_rp)

compare_models_rmse <- 
    model_folds %>% 
    left_join(rmse_estimates) 

rsq_anova <-
    perf_mod(
        compare_models_rmse,
        prior_intercept = student_t(df = 1),
        chains = 4,
        iter = 5000,
        seed = 2
    )

model_post <- 
    rsq_anova %>% 
    # Take a random sample from the posterior distribution
    # so set the seed again to be reproducible. 
    tidy(seed = 35) %>% 
    as_tibble() 

model_post %>% 
    mutate(model = fct_inorder(model)) %>%
    ggplot(aes(x = posterior)) + 
    geom_histogram(bins = 50, col = "white", fill = "blue", alpha = 0.4) + 
    facet_wrap(~ model, ncol = 1) + 
    labs(x = expression(paste("Posterior for RMSE")))

# drawing the feature importance ------------------------------------------
xgb_model <- xgb_model %>% 
    fit(Sale_Price~., data = ames_df)
final_fit <- fit(XGBoost_model, unc13a_cryptic_leaf_psi_log10 ~ ., model_test)
final_rs <- last_fit(XGBoost_model, turbine_capacity ~ ., wind_split)
test_ex <- bake(prep(simple_rec), new_data = model_test)
test_ex$sample = NULL
x_data <- test_ex %>% dplyr::select(-unc13a_cryptic_leaf_psi_log10) %>% as.data.frame()

xgb_predictor <- Predictor$new(XGBoost_res, data = x_data, y = test_ex$unc13a_cryptic_leaf_psi_log10)

xgb_feature_imp <- FeatureImp$new(xgb_predictor, loss = "mae")

rf_feature_imp %>% plot()

xgb_feature_imp %>% plot()

rf_ale <- FeatureEffects$new(rf_predictor)
xgb_ale <- FeatureEffects$new(xgb_predictor)

rf_ale %>% plot()

xgb_ale %>% plot()
