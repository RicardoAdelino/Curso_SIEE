my_pred <- function (mod, envs, pred.type = "cloglog", doClamp = TRUE, ...) {
    requireNamespace("maxnet", quitely = TRUE)
    envs.pts <- as.data.frame(terra::values(envs))
    mxnet.p <- predicts::predict(mod, envs.pts, type = pred.type, clamp = doClamp, 
        ...)
        if(length(dim(mxnet.p)) > 1){
            envs.pts[as.numeric(row.names(mxnet.p)), "pred"] <- mxnet.p
            pred <- terra::rast(cbind(terra::crds(envs, na.rm = FALSE), 
            envs.pts$pred), type = "xyz")
            names(pred) <- "pred"
        return(pred)
        }else {
            mxnet.p <- as.matrix(mxnet.p)
            row.names(mxnet.p) <- 1:length(mxnet.p)
            envs.pts[as.numeric(row.names(mxnet.p)), "pred"] <- mxnet.p
            pred <- terra::rast(cbind(terra::crds(envs, na.rm = FALSE), 
            envs.pts$pred), type = "xyz")
            names(pred) <- "pred"
        return(pred)            
    }
}
