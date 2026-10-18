#devtools::install_github("marlonecobos/kuenm")
#devtools::install_github("johnbaums/rmaxent")
#devtools::install_github("jamiemkass/ENMeval")
#install.packages(
#    c(
#        'devtools',
#        'RColorBrewer', 
#        'geodata', 
#        'rJava', 
#        'ENMeval', 
#        'usdm', 
#        'terra', 
#        'sf', 
#        'janitor',
#        'tidyverse'
#    )
#)

# Carrega os pacotes
gc()
library(tidyverse)
library(sf)
library(terra)
library(sdm)
library(parallel)
library(usdm) 
library(rmaxent)
library(ENMeval)
library(rJava)
library(kuenm)
library(ggplot2)
library(geobr)
library(RColorBrewer)

# Carrega registros de ocorrencia
ocor_ <- readRDS("simposio_db/oco_.rds")

# Prepara função pAUC
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
}

# Carrega dados climáticos
wrl_cl <- geodata::worldclim_global('worldclim', var='bio', res = 5) %>% 
    tidyterra::rename_with(
        ~ str_replace(., "wc2.1_5m_", "")
)

# Extraindo dados
env_data <- terra::extract(
        wrl_cl,
        ocor_[,c("Longitude", "Latitude")],
        xy = TRUE
    ) %>%  
    # remove coluna de ID
    dplyr::select(-ID) %>% 
    # Reajusta posição das coordenadas
    dplyr::relocate(
        tail(names(.), 2), 
        #.before = "bio_1"
        .before = 1
    ) %>% 
    # Remove NA
    tidyr::drop_na() %>% 
    # Renomeia as colunas
    dplyr::rename(
        "long_dec" = "x", 
        "lat_dec" = "y"
)

# Seleção de variaveis por VIF <= 10
var_sel <- usdm::vifstep(
    env_data %>% 
    dplyr::select(
        -long_dec,
        -lat_dec
    ), 
    th = 10
    ) %>%
    purrr::pluck(
        ., 
        "results"
    ) %>% 
dplyr::pull(Variables)

# Direção de fluxo para parametro de raio ao redor da ocorrência
# ==> Na presença de buffer
buf_ <- sf::st_as_sf(
    ocor_ %>% na.omit(),
    coords = c("long_dec", "lat_dec"), 
    remove = FALSE, 
    crs = 4326 
    ) %>% 
    sf::st_transform("ESRI:54009") %>% 
    sf::st_buffer(., dist = 500000) %>% 
sf::st_transform(4326)

# recorta para area de interesse
env <- terra::vect(buf_) %>% terra::crop(wrl_cl, .)

# mascara para area de interesse
env_ <- terra::mask(env, terra::vect(buf_)) %>% terra::subset(.,var_sel)
# ==> Na ausencia de buffer
env_ <- wrl_cl %>% terra::subset(.,var_sel)  

# Background
set.seed(123456)
bg <- env_ %>% 
    terra::spatSample(
        ., 
        size = nrow(ocor_), 
        na.rm = TRUE, 
        values = TRUE, 
        xy = TRUE) %>% 
dplyr::rename("long_dec" = "x", "lat_dec" = "y")

# Prepara os dados pra rodar SDM
occ_data <- env_data %>% 
    dplyr::select(long_dec,lat_dec,var_sel) 

bg_data <- bg 
head(occ_data)

# Run SDM
models_ <- ENMeval::ENMevaluate(
    occs = occ_data[,1:2], 
    bg = bg_data[,1:2],
    envs = env_,
    algorithm =  "maxent.jar", 
    partitions = "block", 
    raster.preds = TRUE,
    #parallel = do_parallel, 
    #numCores = n_cores,
    updateProgress = T, 
    tune.args = list(
        fc =  c("L", "H", "LQ", "LHP", "LQH", "LQHP"),
        rm = seq(.5,10,.5)
    ), 
    user.eval = proc
)

# Visualiza modelos
res <- eval.results(models_)

# Seleciona modelos
opt_aicc <- res %>% dplyr::filter(delta.AICc == 0) %>% 
        dplyr::slice(which.max(auc.val.avg))

# Seleciona parametro
all_pred_ <- models_@models[["fc.H_rm.3"]] 

# Contribuição das variáveis
all_pred_contrib <- all_pred_@results %>% 
    janitor::clean_names() %>%  
    as.data.frame() %>% 
    rownames_to_column(., "betas") %>% 
    as_tibble() %>% 
    filter(
        str_detect(betas, "contribution")
)

# Limiares de corte
thr_all_pred_contrib <- all_pred_@results %>% 
    janitor::clean_names() %>%  
    as.data.frame() %>% 
    rownames_to_column(., "betas") %>% 
    as_tibble() %>% 
    filter(
        str_detect(betas, "maximum_training_sensitivity_plus_specificity_cloglog_threshold")
)

pred.seq <- ENMeval::eval.predictions(res)
plot(pred.seq, col = terrain.colors(100))
