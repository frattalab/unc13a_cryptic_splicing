library(data.table)
library(tidyverse)
library(tidymodels)

library(broom)
library(iml)

clean_data_table = fread(file.path(here::here(),"data","nygc_junction_information.csv"))
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

ranger_recipe <- 
    recipe(formula = unc13a_cryptic_leaf_psi ~ ., data = model_df) %>% 
    step_string2factor(one_of("disease_group", "tissue_clean", "reported_mutations", 
                              "sex")) 

test_ex = bake(prep(ranger_recipe),new_data = model_df)

test_ex = bake(prep(ranger_recipe),new_data = model_df)
linear_model <- lm(unc13a_cryptic_leaf_psi~., data = test_ex)
test_ex %>% 
    cbind(.,.pred = predict(linear_model,test_ex)) %>% 
    ggplot(aes(x = unc13a_cryptic_leaf_psi, y = .pred)) + 
    geom_point() + 
    stat_cor()
#Coef plot
linear_model %>% 
    tidy(conf.int = TRUE) %>% 
    dplyr::slice(-1) %>% 
    mutate(significant_coeff = p.value < 0.05) %>% 
    mutate(term = fct_reorder(term,estimate)) %>% 
    ggplot(aes(x = term, y = estimate, fill = significant_coeff)) + 
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) + 
    geom_point(pch = 21) + 
    coord_flip() 



# Start making some non-interpretable models 

random_forest_model <- rand_forest() %>%
    set_mode("regression") %>%
    set_engine("randomForest")

xgb_model <- boost_tree() %>%
    set_mode("regression") %>%
    set_engine("xgboost")




#Train models

set.seed(42)


random_forest_model <- rand_forest() %>% 
    set_mode("regression") %>% 
    set_engine("ranger") %>% 
    fit(unc13a_cryptic_leaf_psi~. , data = test_ex)

test_ex %>% 
    cbind(.,predict(random_forest_model,test_ex)) %>% 
    ggplot(aes(x = unc13a_cryptic_leaf_psi, y = .pred)) + 
    geom_point() + 
    stat_cor() + ggtitle("An OverFit Random Forest")

xgb_model <- boost_tree() %>% 
    set_mode("regression") %>% 
    set_engine("xgboost") %>% 
    fit(unc13a_cryptic_leaf_psi~., data = test_ex) 

test_ex %>% 
    cbind(.,predict(xgb_model,test_ex)) %>% 
    ggplot(aes(x = unc13a_cryptic_leaf_psi, y = .pred)) + 
    geom_point() + 
    stat_cor() + ggtitle("An OverFit XGBoost")

x_data <- test_ex %>% dplyr::select(-unc13a_cryptic_leaf_psi) %>% as.data.frame()

xgb_predictor <- Predictor$new(xgb_model, data = x_data, y = test_ex$unc13a_cryptic_leaf_psi)

rf_predictor <- Predictor$new(random_forest_model, data = x_data, y = test_ex$unc13a_cryptic_leaf_psi)

rf_feature_imp <- FeatureImp$new(rf_predictor, loss = "rmse")

xgb_feature_imp <- FeatureImp$new(xgb_predictor, loss = "rmse")

rf_feature_imp %>% plot() + ggtitle("Random Forest Model - Change in RMSE") + xlim(1,2)

xgb_feature_imp %>% plot() + ggtitle("XGBoost Model - Change in RMSE") + xlim(1,2.5)

rf_ale <- FeatureEffects$new(rf_predictor)
xgb_ale <- FeatureEffects$new(xgb_predictor)

rf_ale %>% plot(ncol = 2)

xgb_ale %>% plot(ncol = 2)
























