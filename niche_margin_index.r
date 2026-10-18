# Load packages
gc()
library(tidyverse)
library(terra)
library(ade4)
library(ecospat)
library(ks)
library(sp)
library(sf)
library(usdm)

# Generate climatic data at the chosen grain
clim <- geodata::worldclim_global('worldclim', var='bio', res = 5)

#Load South America Non-Native occurrences
inv_occ <- readRDS("simposio_db/oco_.rds") %>%
    dplyr::filter(inAF == 0) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(
        sci_name, 
        Longitude, 
        Latitude, 
        name_long
)

# Keep Country flag
inv_occ_env <- terra::extract(
        clim,
        inv_occ[,c("Longitude","Latitude")]
    ) %>% 
    dplyr::select(-1) %>% 
tibble::add_column(inv_occ[,4])

#Load Native range
all_native <- readRDS("simposio_db/oco_.rds") %>%
    dplyr::filter(inAF == 1) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(
        sci_name, 
        Longitude, 
        Latitude, 
        name_long
)

nat_occ_env <- terra::extract(
        clim,
        all_native[,c("Longitude","Latitude")]
    ) %>% 
    dplyr::select(-1) 
#tibble::add_column(all_native[,4]) 

# NMI routine
# Select variables
sp_vars <- dplyr::bind_rows(nat_occ_env,inv_occ_env) %>% 
    dplyr::select(-last_col()) %>% 
usdm::vifstep(., th =10) 

# Keep var vectors    
sp_vars_name <- sp_vars@results$Variables

# Create bakround variables
clim_df <- clim %>% 
    tidyterra::select(sp_vars_name) %>% 
    base::as.data.frame() %>% 
stats::na.exclude()

# Create native subset
nat_sub <- nat_occ_env %>% select(sp_vars_name)

# Create invaded subset
inv_sub <- inv_occ_env %>% select(sp_vars_name)

# Run Component Analisys in Bkacground data
pca <- ade4::dudi.pca(clim_df,scannf = FALSE, center = FALSE, nf = 2)

# Keep Background scores
bkg_scores <- pca$li

# Smooth the bakground densities
fhat <- ks::kde(bkg_scores, compute.cont = TRUE)

# Keep contour lines of background densities (sometimes fail)
c99 <- contourLines(
    fhat$eval.points[[1]],
    fhat$eval.points[[2]],
    fhat$estimate,
    level = 0.0001
    #level = quantile(fhat$estimate, 0.9, na.rm = TRUE)
)

# Create the polygon of contour
l99 <- list()        
for(k in 1:length(c99)){
   l99[[k]] = Polygons(list(Polygon(c99[[k]][-1])),ID=k)
}

bkg_pol <- SpatialPolygons(l99)
bkg_pol <- aggregate(bkg_pol)

# species scores
sp_scores <- na.omit(suprow(pca,nat_sub)$li)

fhat <- kde(sp_scores,compute.cont=TRUE)
c99 <- contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=fhat$cont[99])

l99<-list()        
for(k in 1:length(c99)){
    l99[[k]] <- Polygons(list(Polygon(c99[[k]][-1])),ID=k)  
}

sp_pol <- SpatialPolygons(l99)
sp_pol <- raster::aggregate(sp_pol)

# invasion scores
intros_scores <- suprow(pca,inv_sub)$li
NArows <- is.na(intros_scores$Axis1)
intros <- inv_sub[!NArows,]
intros_scores <- intros_scores[!NArows,]
coordinates(intros_scores) <- cbind(intros_scores$Axis1 , intros_scores$Axis2) # DataFrame to SpatialPointsDataFrame

### calculate NMI values
NMI <- function(foc.pop,niche){
  margin<-as(niche, 'SpatialLines') #niche margin
  niche.regpoints <- spsample(niche,n=100000,"regular") #regular points within the niche
  max.inner <- max(sf::st_distance(st_as_sf(niche.regpoints),st_as_sf(margin)))
  
  inout<-over(foc.pop,niche) # belonging to the niche
  NMI.abs <- sf::st_distance(st_as_sf(foc.pop), st_as_sf(niche))
  #NMI.abs<-gDistance(spgeom1=foc.pop, spgeom2=margin, byid=T)
  NMI.abs[which(is.na(inout))]<- -NMI.abs[which(is.na(inout))]
  NMI<-NMI.abs/max.inner
  
  return (list(NMI = c(NMI),NMI.abs=c(NMI.abs)))
}

bkg_grid <- raster::raster(nrows=1000,ncols=1000,raster::extent(bkg_pol))
bkg <- raster::rasterize(bkg_pol,bkg_grid)
bkg_pts <- data.frame(coordinates(bkg))
coordinates(bkg_pts) <- cbind(bkg_pts$x , bkg_pts$y)

bkg_NMI <- NMI(foc.pop = bkg_pts,niche = sp_pol)
intro_NMI <- NMI(foc.pop = intros_scores,niche=sp_pol)

#plot NMI
in_bkg <- !is.na(values(bkg)) #pixels in bkg
pos <- bkg_NMI$NMI>=0 #pixels with positive values

#plot innerness
inner <- bkg  
inner[in_bkg&pos] <- bkg_NMI$NMI[in_bkg&pos] 
inner[inner ==1] <- NA
Pal_inner <- colorRampPalette(c('white', 'steelblue', 'blue'))

#plot outerness
outer <- bkg  
outer[in_bkg&!pos] <- bkg_NMI$NMI[in_bkg&!pos]
outer[outer == 1] <- NA
Pal_outer <- colorRampPalette((c('white','orange')))

# point
inv_occ_country <- inv_occ_env %>% pull(name_long)

point_label <- intros_scores@data %>% 
    add_column(
        country = inv_occ_country,
        NMI = intro_NMI[[1]]
    ) %>% 
drop_na(country)    

point_label_ <- point_label %>% distinct(Axis1,Axis2,.keep_all = TRUE)

{x11()}
# plot axes
plot(coordinates(bkg_pts),type="n",xlab="PC1",ylab="PC2", xlim = c(-25,15))
plot(inner,col=Pal_inner(100),add=T, legend=FALSE)
plot(outer,col=Pal_outer(100),add=T, legend=FALSE)
# accessible area
plot(bkg_pol,add=T,lty=1) 
# native climatic niche margin
plot(sp_pol,add=T,border = "red") 


points(
    point_label_ %>% filter(country == "Uruguay"),
    pch = 19,
    col = 'red'
)

points(
    point_label_ %>% filter(country == "Uruguay"),
    pch = 21
)
