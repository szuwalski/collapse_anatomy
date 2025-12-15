#===============FISHERIES data=========================
#==plot maps of effort and CPUE
#==plot discards and catches?



#===============SURVEY DATA==========================
#==plot map of legal and large male centroids and distributions
#==plot time series
#==plot management advice by different bodies
library(dplyr)
library(reshape2)
library(ggplot2)
library(png)
library(grid)
library(patchwork)
library(crabpack)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggridges)
library(tidyr)
library(rnaturalearth)

## Pull specimen data
specimen_data <- crabpack::get_specimen_data(species = "SNOW",
                                             region = "EBS",
                                             years = c(1982:2025),
                                             channel = 'API')

male_snow_ind_101 <- crabpack::calc_bioabund(crab_data = specimen_data,
                                             species = "SNOW",
                                             region = "EBS",
                                             sex='male',
                                             size_min = 101)
male_snow_ind_101$size<-">101mm"
male_snow_ind_all <- crabpack::calc_bioabund(crab_data = specimen_data,
                                             species = "SNOW",
                                             region = "EBS",
                                             sex='male',
                                             size_min=40)
male_snow_ind_all$size<-">40mm"
male_snow_spatial_101 <- crabpack::calc_cpue(crab_data = specimen_data,
                                             species = "SNOW",
                                             region = "EBS",
                                             sex='male',
                                             size_min = 101)

male_snow_ind_76 <- crabpack::calc_bioabund(crab_data = specimen_data,
                                            species = "SNOW",
                                            region = "EBS",
                                            sex='male',
                                            size_min = 76)

male_snow_spatial_76 <- crabpack::calc_cpue(crab_data = specimen_data,
                                            species = "SNOW",
                                            region = "EBS",
                                            sex='male',
                                            size_min = 76)

male_snow_all_size<-crabpack::calc_bioabund(crab_data = specimen_data,
                                            species = "SNOW",
                                            region = "EBS",
                                            sex='male',
                                            size_min=40,bin_1mm=TRUE)

dat_in<-rbind(male_snow_ind_101,male_snow_ind_all)
dat_in$upci<-dat_in$ABUNDANCE * exp(1.96*sqrt(log(1+dat_in$ABUNDANCE_CV^2)))
dat_in$dnci<-dat_in$ABUNDANCE / exp(1.96*sqrt(log(1+dat_in$ABUNDANCE_CV^2)))

traj_plot<-ggplot()+
  geom_point(data=filter(dat_in,YEAR<2026),aes(x=YEAR,y=ABUNDANCE))+
  geom_errorbar(data=filter(dat_in,YEAR<2026),aes(x=YEAR,ymin=dnci,ymax=upci))+
  facet_wrap(~size,ncol=1,scales='free_y')+
  theme(legend.position='none')+theme_bw()+
  expand_limits(y=0)+xlim(1980,2026)

take_yr<-seq(2016,2026)
for(x in 1:length(take_yr))
{
  if(take_yr[x]!=2020)
  {
  traj_plot<-ggplot()+
  geom_point(data=filter(dat_in,YEAR<take_yr[x]),aes(x=YEAR,y=ABUNDANCE))+
    geom_line(data=filter(dat_in,YEAR<take_yr[x]),aes(x=YEAR,y=ABUNDANCE),col='blue',lwd=1.2)+
  geom_errorbar(data=filter(dat_in,YEAR<take_yr[x]),aes(x=YEAR,ymin=dnci,ymax=upci))+
  facet_wrap(~size,ncol=1,scales='free_y')+
  theme(legend.position='none')+theme_bw()+
  expand_limits(y=0)+xlim(1980,2026)

  indat<-melt(select(male_snow_all_size,c("YEAR","SIZE_1MM","ABUNDANCE")),id.vars=c("YEAR","SIZE_1MM"))
  p <- ggplot(dat=filter(indat,YEAR<take_yr[x])) 
  p <- p + geom_density_ridges(aes(x=SIZE_1MM, y=YEAR, height = value,
                                   group = YEAR, 
                                   alpha=.9999),fill='blue',stat = "identity",scale=5) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90)) +
    labs(x="Carapace width (mm)",y="Year") +
    xlim(40,130)+ylim(1980,2027)
  
  layouts<-"
    122
    "
  outs<-p+traj_plot+plot_layout(design=layouts)
  
  
  png(paste("plots/traj",take_yr[x],".png"),height=8,width=12,res=350,units='in') 
  print(outs)
  dev.off()
  }
}
###############################################################
## 1. Libraries
###############################################################

library(sf)
library(dplyr)
library(purrr)
library(ggplot2)
library(ks)


###############################################################
## 2. Load data
## df must have LATITUDE, LONGITUDE, CPUE_MT, YEAR
###############################################################

df <- male_snow_spatial_76   # <-- replace with your object


###############################################################
## 3. Convert to sf and project to EPSG 3338 (Alaska Albers)
###############################################################

df_sf <- st_as_sf(df,
                  coords = c("LONGITUDE","LATITUDE"),
                  crs = 4326)

df_proj <- st_transform(df_sf, 3338)


###############################################################
## 4. KDE → 50% and 90% UD contours using ks::kde
###############################################################

compute_kde_contours <- function(sf_year,
                                 gridsize = c(120, 120)) {
  
  coords <- st_coordinates(sf_year)
  w      <- sf_year$CPUE_MT
  
  ok <- is.finite(coords[,1]) & is.finite(coords[,2]) & is.finite(w)
  coords <- coords[ok,,drop=FALSE]
  w      <- w[ok]
  
  if (nrow(coords) < 5 || sum(w) == 0)
    return(NULL)
  
  ## KDE (ks chooses bandwidth unless you provide H)
  kde_res <- kde(
    x        = coords,
    w        = w,
    gridsize = gridsize
  )
  
  xg <- kde_res$eval.points[[1]]
  yg <- kde_res$eval.points[[2]]
  z  <- kde_res$estimate
  
  ## Normalize to probability surface
  p <- z / sum(z, na.rm=TRUE)
  
  ## UD thresholds
  v <- as.vector(p)
  o <- sort(v, decreasing = TRUE)
  cumu <- cumsum(o)
  
  thr50 <- o[which(cumu >= 0.25)[1]]
  thr90 <- o[which(cumu >= 0.50)[1]]
  
  cls <- contourLines(x = xg, y = yg, z = p,
                      levels = c(thr50, thr90))
  
  if (length(cls) == 0) return(NULL)
  
  polys <- lapply(cls, function(cl) {
    pts <- cbind(cl$x, cl$y)
    if (!all(pts[1,] == pts[nrow(pts),]))
      pts <- rbind(pts, pts[1,])
    st_polygon(list(pts))
  })
  
  levels <- sapply(cls, `[[`, "level")
  pct    <- ifelse(abs(levels - thr50) < 1e-12, 25, 50)
  
  sf_obj <- st_sf(
    percent  = pct,
    geometry = st_sfc(polys, crs = st_crs(sf_year))
  ) %>%
    group_by(percent) %>%
    summarise(geometry = st_union(geometry), .groups="drop")
  
  sf_obj
}


###############################################################
## 5. Compute contours for each year
###############################################################

contours <- df_proj %>%
  split(.$YEAR) %>%
  imap_dfr(function(sf_year, yr) {
    out <- compute_kde_contours(sf_year)
    if (!is.null(out) && nrow(out) > 0) {
      out$YEAR <- yr
      out
    }
  })


###############################################################
## 6. Compute biomass-weighted centroids
###############################################################

coords_all <- st_coordinates(df_proj)

centroids <- df_proj %>%
  mutate(x = coords_all[,1],
         y = coords_all[,2]) %>%
  group_by(YEAR) %>%
  summarise(
    cx = weighted.mean(x, CPUE_MT, na.rm=TRUE),
    cy = weighted.mean(y, CPUE_MT, na.rm=TRUE)
  ) %>%
  st_as_sf(coords = c("cx","cy"), crs = 3338)


###############################################################
## 7. Anti-meridian-safe transformation to lon/lat
###############################################################

fix_antimeridian <- function(sf_obj) {
  
  # Step 1: Use a lon_wrap CRS (keeps polygons continuous)
  ll360 <- "+proj=longlat +datum=WGS84 +lon_wrap=180"
  obj360 <- st_transform(sf_obj, ll360)
  
  # Step 2: Split across the 180° line
  obj_split <- st_wrap_dateline(obj360,
                                options = c("WRAPDATELINE=YES",
                                            "DATELINEOFFSET=180"),
                                quiet = TRUE)
  
  # Step 3: Transform back to standard WGS84 lon/lat
  st_transform(obj_split, 4326)
}

contours_ll  <- fix_antimeridian(contours)
centroids_ll <- st_transform(centroids, 4326)


###############################################################
## 8. Optional: clip contours to buffered hull of all points
## (useful to prevent “tails” into low density areas)
###############################################################

hull <- df_proj %>%
  st_union() %>%
  st_convex_hull() %>%
  st_buffer(75000)   # 75 km buffer

hull_ll <- st_transform(hull, 4326)

contours_ll <- st_intersection(contours_ll, hull_ll)


###############################################################
## 9. Plot (no points)
###############################################################
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<- -179
lon_2<- -156
lat_1<- 52
lat_2<- 65

png(paste("plots/densities_101.png"),height=10,width=8,res=350,units='in')
ggplot() +
  geom_sf(data = filter(contours_ll,YEAR>2016),
          aes(color = as.factor(YEAR),
              lty= as.factor(percent)),fill=NA,  alpha = 0.25, linewidth=1) +
 # geom_sf(data = centroids_ll,
 #         aes(color = as.factor(YEAR)),
 #         size = 3) +
 # scale_fill_viridis_d(option="C", name = "Year") +
  scale_color_viridis_d(option="D", name = "Percent / Year") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  labs(title = "Eastern Bering Sea Snow Crab",
       subtitle = "Biomass-weighted centroids and 25%/50% KDE utilization distributions",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none")+facet_wrap(~YEAR)
dev.off()
