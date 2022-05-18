# ---------------------------------------
# Brainplot functions!
#
# These take as arguments:
#  dta: your dataframe of values and what regions they're assigned to
#  out_prefix: what should be the start of your output file name
#  out_suffix: what goes at the end of your output file name (after the region and looping variable info)
#  loop_over = the name of the column that contains a variable that you want to loop over
#     (e.g. make brainplots for each of several p-value thresholds). Can be a dummy variable with only 1 value.
#  region_col = what column has the brain regions for the data, as named in the "regionordering" variable below
#  Z_col = the column with the numeric data you want to plot
#  analysis_col = This will go into the title of the plots, contents should be things like "Surface Area" or "Thickness"
#  max_val = the highest value in your Z_col, you may want to make this consistent across several plots
#  low_color = bottom of the color scale
#  high_color = top of the color scale
#
# The output of these functions is an HTML file that you can open in a web browswer (it will likely crash Firefox, use Chrome).
# From the open browser window, there's a menu on the top right that includes a save to PNG option. This will give you files
# ready to go into your figures.
#
# One day if I have time I will expand these a little and put them into a nice package.
# Note: credit for the original versions goes to Jason Stein, PhD at Univ. of North Carolina Chapel Hill.
# -----------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(here)
library(GenomicRanges)
library(biomaRt)
library(diffloop)
library(plotly)


### Load descriptions of the brain surface and regions
#regionordering <- read.csv(here("plotting", "freesurfer_orderandcolor.csv"))
#load(here("plotting", "FreesurferRegionalObjs.Rdata"))
regionordering <- read.csv("/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/Cerebral_Cortex_revisions/plotting/freesurfer_orderandcolor.csv")
load("/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/Cerebral_Cortex_revisions/plotting/FreesurferRegionalObjs.Rdata")
str(objs)




### Surface Area ###
# This function is fully commented, the Thickness version is essentially the same
brainplot_SA <- function(dta, 
                         out_prefix, 
                         out_suffix, 
                         loop_over, 
                         region_col, 
                         Z_col, 
                         analysis_col, 
                         max_val, 
                         low_color, 
                         high_color, 
                         nonsig_color) {

  ## Load your data frame and set columns to be used for plotting
  fZ <- dta
  fZ$SNP <- fZ[, loop_over]
  fZ$Trait <- fZ[, region_col]

  fZ$Z <- fZ[, Z_col]
  Z.SA <- fZ[grepl("Surf", fZ[[analysis_col]]), ] # sets a smaller dataframe of just the surface area data
  Z.TH <- fZ[grepl("Thick", fZ[[analysis_col]]), ] # does the same, but for thickness data

  ## Cleaning up the names of the brain regions in the objs list, and removing the left hemisphere versions (we only need one set)
  objs <- objs[1:35]
  names(objs) <- sapply(names(objs), function(x) {
    unlist(strsplit(x, ".", fixed = TRUE))[length(unlist(strsplit(x, ".", fixed = TRUE)))]
  })

  ## Settings to remove axes from the plots
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )

  # -----------------------------------
  # make the plots
  # -----------------------------------
  my_colors = colorRampPalette(c(low_color,"white", high_color))(50)
  my_colors[25] = nonsig_color

  # Setting a min/max threshold for the color labeling that will be consistent across all
  # annotations being plotted, since we only subset by annotation within the plotting for-loop
  # Z.SA.max = max(Z.SA$Z)
  # Z.TH.max = max(Z.TH$Z)

  # Setting up a list to loop over, in case the input dataframe has data for multiple plots
  uniquesnps.SA <- unique(Z.SA$SNP)
  uniquesnps.TH <- uniquesnps.SA
  
  ## Looping over the loop variable, making a pair of medial/lateral brain plots each time
  for (j in 1:length(uniquesnps.SA)) {
    thissnp <- uniquesnps.SA[[j]]
    snpind <- which(Z.SA$SNP == thissnp)
    thisSA <- Z.SA[snpind, ]
    # fullSA.padj <- thisSA$fdr[1]
    # fullSA.Z <- thisSA$Z[1]

    ## Make sure the order of objs (the brain surfaces we loaded), matches the order of the input data
    inds <- match(names(objs), thisSA$Trait)
    thisSA <- thisSA[inds, ]

    ## Ensure the color scale is symmetrical from positive to negative
    maxZ <- max_val
    minZ <- -(max_val)
    indmax <- which.max(c(abs(maxZ), abs(minZ)))
    if (indmax == 1) {
      minZ <- -abs(maxZ)
    } else {
      maxZ <- abs(minZ)
    }

    ## Loop through each brain region for plotting, using plotly's "mesh3d" style
    for (i in 1:length(objs)) {
      ftitle <- paste("SA", thissnp) # title of the plot
      obj <- objs[[i]]
      ## Build up the plot, one region at a time (add_mesh adds new regions to an existing plot)
      if (i == 1) {
        p <- plot_ly(
          type = "mesh3d",
          x = obj$shapes[[1]]$positions[1, ], # x coord
          y = obj$shapes[[1]]$positions[2, ], # y coord
          z = obj$shapes[[1]]$positions[3, ], # z coord
          i = obj$shapes[[1]]$indices[1, ], # vertex indices, the first vertex of the triangle
          j = obj$shapes[[1]]$indices[2, ], # vertex indices, the second vertex of the triangle
          k = obj$shapes[[1]]$indices[3, ], # vertex indices, the third vertex of the triangle
          intensity = rep(thisSA$Z[i], length(obj$shapes[[1]]$positions[1, ])), # sets the value to plot on surface of the mesh
          opacity = 1, # opacity of the trace
          cauto = FALSE, # should we set the color range based on the values for this region?
          cmax = maxZ, # what value is the top of the color range?
          cmin = minZ, # bottom of color range
          #colors = colorRampPalette(c(low_color, "grey", high_color))(50), # use 50 breaks in the diverging color palette
          colors = my_colors,
          showscale = TRUE
        ) # show the color scale
      } else {
        p <- add_mesh(
          p = p,
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          intensity = rep(thisSA$Z[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          showscale = FALSE
        )
      }
    }
    # Save an html file where the "camera" is pointed at the medial surface of the brain
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = 2, y = 0, z = 0))), title = ftitle)
    htmlwidgets::saveWidget(p, file = paste0(out_prefix, uniquesnps.SA[j], "_medial_SA_", out_suffix, ".html"), selfcontained = FALSE)
    # Same for the lateral surface of the brain
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = -2, y = 0, z = 0))), title = ftitle)
    htmlwidgets::saveWidget(p, file = paste0(out_prefix, uniquesnps.SA[j], "_lateral_SA_", out_suffix, ".html"), selfcontained = FALSE)
  }
}


### Thickness ###

brainplot_TH <- function(dta, out_prefix, out_suffix, loop_over, region_col, Z_col, analysis_col, max_val, low_color, high_color) {
  fZ <- dta
  fZ$SNP <- fZ[, loop_over]
  fZ$Trait <- fZ[, region_col]
  fZ$Z <- fZ[, Z_col]
  Z.SA <- fZ[grepl("Surf", fZ[[analysis_col]]), ]
  Z.TH <- fZ[grepl("Thick", fZ[[analysis_col]]), ]

  objs <- objs[1:35]
  names(objs) <- sapply(names(objs), function(x) {
    unlist(strsplit(x, ".", fixed = TRUE))[length(unlist(strsplit(x, ".", fixed = TRUE)))]
  })

  ## Remove axes
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )

  # -----------------------------------
  # make the plots
  # -----------------------------------

  # Setting a min/max threshold for the enrichment color labeling that will be consistent across all
  # annotations being plotted, since we only subset by annotation within the plotting for-loop

  uniquesnps.SA <- unique(Z.SA$SNP)
  uniquesnps.TH <- uniquesnps.SA

  for (j in 1:length(uniquesnps.TH)) {
    thissnp <- uniquesnps.TH[[j]]
    snpind <- which(Z.TH$SNP == thissnp)
    thisTH <- Z.TH[snpind, ]
    fullTH.padj <- thisTH$fdr[1]
    fullTH.Z <- thisTH$Z[1]
    ## Make sure the order of each corresponds
    inds <- match(names(objs), thisTH$Trait)
    thisTH <- thisTH[inds, ]
    
    ## Ensure the color scale is symmetrical from positive to negative
    maxZ <- max_val
    minZ <- -(max_val)
    indmax <- which.max(c(abs(maxZ), abs(minZ)))
    if (indmax == 1) {
      minZ <- -abs(maxZ)
    } else {
      maxZ <- abs(minZ)
    }
    ## Loop through each region for plotting
    for (i in 1:length(objs)) {
      ftitle <- paste("TH", thissnp)
      obj <- objs[[i]]
      if (i == 1) {
        p <- plot_ly(
          type = "mesh3d",
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          intensity = rep(thisTH$Z[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          colors = colorRampPalette(c(low_color, "grey", high_color))(50),
          showscale = TRUE
        )
      } else {
        p <- add_mesh(
          p = p,
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          intensity = rep(thisTH$Z[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          showscale = FALSE
        )
      }
    }
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = 2, y = 0, z = 0))), title = ftitle)
    htmlwidgets::saveWidget(p, file = paste0(out_prefix, uniquesnps.TH[j], "_medial_TH_", out_suffix, ".html"), selfcontained = FALSE)
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = -2, y = 0, z = 0))), title = ftitle)
    htmlwidgets::saveWidget(p, file = paste0(out_prefix, uniquesnps.TH[j], "_lateral_TH_", out_suffix, ".html"), selfcontained = FALSE)
  }
}

### Full Surface Area or Thickness ###
brainplot_full <- function(dta, out_prefix, out_suffix, loop_over, region_col, Z_col, analysis_col, max_val, low_color, high_color, nonsig_color) {
  fZ <- dta
  fZ$SNP <- fZ[, loop_over]
  fZ$Trait <- fZ[, region_col]

  fZ$Z <- fZ[, Z_col]
  Z.SA <- fZ[grepl("Surf", fZ[[analysis_col]]), ]
  Z.TH <- fZ[grepl("Thick", fZ[[analysis_col]]), ]

  columns <- c("Trait", "SNP", "Z")
  FullTH <- Z.TH[Z.TH$Region == "Full", columns]
  Z.TH.Full <- left_join(Z.TH, FullTH, by = "SNP")

  FullSA <- Z.SA[Z.SA$Region == "Full", columns]
  Z.SA.Full <- left_join(Z.SA, FullSA, by = "SNP")

  ## Take only the left hemisphere obj files (everything in ENIGMA is mean bilateral)
  ## So no need for both hemispheres
  objs <- objs[1:35]
  names(objs) <- sapply(names(objs), function(x) {
    unlist(strsplit(x, ".", fixed = TRUE))[length(unlist(strsplit(x, ".", fixed = TRUE)))]
  })

  ## Remove axes
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )

  # -----------------------------------
  # make the plots
  # -----------------------------------
  my_colors = colorRampPalette(c(low_color,"white", high_color))(50)
  my_colors[25] = nonsig_color

  # Setting a min/max threshold for the enrichment color labeling that will be consistent across all
  # annotations being plotted, since we only subset by annotation within the plotting for-loop
  # Z.SA.max = max(Z.SA$Z)
  # Z.TH.max = max(Z.TH$Z)

  uniquesnps.SA <- unique(Z.SA$SNP)
  uniquesnps.TH <- uniquesnps.SA

  for (j in 1:length(uniquesnps.TH)) {
    thissnp <- uniquesnps.TH[[j]]
    snpind <- which(Z.TH.Full$SNP == thissnp)
    thisTH <- Z.TH.Full[snpind, ]
    fullTH.padj <- thisTH$fdr[1]
    fullTH.Z <- thisTH$Z[1]
    ## Make sure the order of each corresponds
    inds <- match(names(objs), thisTH$Trait.x)
    thisTH <- thisTH[inds, ]
    maxZ <- max_val
    minZ <- -(max_val)
    ## Ensure the color scale is symmetrical from positive to negative
    indmax <- which.max(c(abs(maxZ), abs(minZ)))
    if (indmax == 1) {
      minZ <- -abs(maxZ)
    } else {
      maxZ <- abs(minZ)
    }
    ## Loop through each region for plotting
    for (i in 1:length(objs)) {
      ftitle <- paste("TH", thissnp)
      obj <- objs[[i]]
      ## If it is the first plot, set up the mesh
      if (i == 1) {
        p <- plot_ly(
          type = "mesh3d",
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          ## for two of the vertices, do the whole range of colors (you can't see the in the plots, but without it I couldn't get the color ranges to work)
          intensity = rep(thisTH$Z.y[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          colors = my_colors,
          showscale = TRUE
        )
      } else {
        p <- add_mesh(
          p = p,
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          intensity = rep(thisTH$Z.y[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          showscale = FALSE
        )
      }
    }
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = -2, y = 0, z = 0))), title = ftitle)
    p_file <- paste0(out_prefix, uniquesnps.TH[j], "_lateral_fullTH_", out_suffix, ".html")
    message(p_file)
    htmlwidgets::saveWidget(p, file = p_file, selfcontained = FALSE)
  }

  for (j in 1:length(uniquesnps.SA)) {
    thissnp <- uniquesnps.SA[[j]]
    snpind <- which(Z.SA.Full$SNP == thissnp)
    thisSA <- Z.SA.Full[snpind, ]
    fullSA.padj <- thisSA$fdr[1]
    fullSA.Z <- thisSA$Z[1]
    ## Make sure the order of each corresponds
    inds <- match(names(objs), thisSA$Trait.x)
    thisSA <- thisSA[inds, ]
    maxZ <- max_val
    minZ <- -(max_val)
    ## Ensure the color scale is symmetrical from positive to negative
    indmax <- which.max(c(abs(maxZ), abs(minZ)))
    if (indmax == 1) {
      minZ <- -abs(maxZ)
    } else {
      maxZ <- abs(minZ)
    }
    ## Loop through each region for plotting
    for (i in 1:length(objs)) {
      ftitle <- paste("SA", thissnp)
      obj <- objs[[i]]
      ## If it is the first plot, set up the mesh
      if (i == 1) {
        p <- plot_ly(
          type = "mesh3d",
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          ## for two of the vertices, do the whole range of colors (you can't see the in the plots, but without it I couldn't get the color ranges to work)
          intensity = rep(thisSA$Z.y[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          colors = my_colors,
          showscale = TRUE
        )
      } else {
        p <- add_mesh(
          p = p,
          x = obj$shapes[[1]]$positions[1, ],
          y = obj$shapes[[1]]$positions[2, ],
          z = obj$shapes[[1]]$positions[3, ],
          i = obj$shapes[[1]]$indices[1, ],
          j = obj$shapes[[1]]$indices[2, ],
          k = obj$shapes[[1]]$indices[3, ],
          intensity = rep(thisSA$Z.y[i], length(obj$shapes[[1]]$positions[1, ])),
          opacity = 1,
          cauto = FALSE,
          cmax = maxZ,
          cmin = minZ,
          showscale = FALSE
        )
      }
    }
    p <- layout(p = p, scene = list(xaxis = ax, yaxis = ax, zaxis = ax, camera = list(eye = list(x = -2, y = 0, z = 0))), title = ftitle)
    htmlwidgets::saveWidget(p, file = paste0(out_prefix, uniquesnps.SA[j], "_lateral_fullSA_", out_suffix, ".html"), selfcontained = FALSE)
  }
}