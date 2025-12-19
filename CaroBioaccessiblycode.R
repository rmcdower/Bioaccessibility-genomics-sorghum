##install packages

install.packages("lme4")
install.packages("qtl")
install.packages("devtools")
library(devtools)
install.packages("qtlDesign", repos="http://R-Forge.R-project.org")
install.packages("remotes")
remotes::install_github("jtlovell/qtlTools")

##load libraries
library(lme4)
library(qtl)
library(qtlDesign)
library(qtlTools)

#load cross
setwd("~/Desktop/qtl")
mydata <- read.cross("csv", file= "AllPhenosQTL.csv", genotypes= c("A","H", "B"), alleles= c("A", "B"))

#convert to RILS
mydata<-convert2riself(mydata)
summary(mydata)
mydata<-est.rf(mydata)
summary(mydata)
plot(mydata)


error<-calc.errorlod(mydata, error.prob = 0.01)

#wont plot all individuals (ind) at same time, look through and check for recombination errors
plotGeno(error, chr= 3, ind= 1:40)
plotGeno(error, ind= 1:40)

#QC
geno.image(mydata)
plotRF(mydata)
errors<-top.errorlod(error, cutoff= 4)
errors

qtlfile <- calc.genoprob(mydata, step=0, error.prob = 0.001)
phenames(qtlfile)
phes <- phenames(qtlfile)

# Then regenerate permutations and run calcCis
modelperms <- scanone(qtlfile, pheno.col=1:31, method="hk", n.perm=1000)
s1<-scanone(qtlfile, pheno.col = 1:31, method = "hk")
cis <- calcCis(qtlfile, s1.output=s1, perm.output=modelperms, lodcolumn=1)

cis

#recheck for normal dist
hist(mydata$pheno$ACrypto22)

#phenotypes columns
phes

#calc credible intervals
s2<-scanone(qtlfile, pheno.col = "TKW23", method = "hk")
#bayes int (scanone object, chr, significance)
bayesout<-bayesint(s2, 3, 0.95, expandtomarkers=TRUE)
plot(bayesout)
bayesout




#look for a marker and make effect plot of alleles
max(s1$lod)
find.marker(mydata,10, 16)
effectplot(mydata,pheno.col = "betaalpha", mname1 = "S10_1493293")
effectplot(mydata,pheno.col = "TKW23", mname1 = "S3_4244979")

phes
summary(modelperms, alpha=c(0.01, 0.05))

summary(s1, perms=modelperms, alpha=0.05, pvalues=TRUE)

pullSigQTL(qtlfile, pheno.col=2 ,s1.output = s1,  perm.output = modelperms)



mods <- pullSigQTL(qtlfile, s1.output = s1,  perm.output = modelperms)
print(mods)
mods<-mods[sapply(mods,length)!=1] # a NULL QTL model has length = 1, toss these
print(mods)

#We can also plot confidence intervals directly from scanone data,
#however, it is important to note that multiple peaks on a single chromosome can violate assumptions of such intervals.
cis<-calcCis(qtlfile,s1.output=s1, perm.output=modelperms)
head(cis)


cis

#adjust colors based on amount of traits mapping
col= c("red", "blue", "green", "purple","darkgreen")
col=c("blue", "purple", "darkgreen", "pink")
col="blue"
segmentsOnMap(qtlfile, calcCisResults = cis, legendCex = .6, palette= rainbow,lwd = 3,chrBuffer = c(.2,.2))


stats<-lapply(names(mods), function(x){
  mod<-mods[[x]]
  form<-paste("y ~",paste(mod$altname, collapse = " + "))
  mod<-refineqtl(qtlfile, qtl = mod, formula = form, 
                 verbose=F, method="hk", pheno.col = 1:31)
  qtlStats(qtlfile, 
           pheno.col = 1:31, 
           form = form, 
           mod = mod)
})
stats <- do.call(rbind, stats)
print(head(stats))

write.table(stats, "statTKW23.csv", sep = ",", row.names = FALSE)

segmentsOnMap(cross=qtlfile, phe=phes,
              col = col, ##add more colors if more phenotypes
              chr=stats[[1]]$chr, 
              l = stats[[1]]$lowposition, 
              lwd=5,
              chrBuffer = c(0.2, 0.2),
              h = stats[[1]]$highposition,
              legendPosition = "bottomright",
              legendCex= 1,
              leg.lwd = 5,
              leg.inset=.05)

with(cis, segmentsOnMap(qtlfile, phe = phes, chr = chr,
                        l = lowposition, h = highposition, legendCex = .5,
                        peakcM = pos,
                        tick.width = .1,  chrBuffer = c(.15,.2), palette = rainbow))
library(paletteer)
library(RColorBrewer)

palette
col= c("coral" , "orange", "orange" , "red", "blue" , "aquamarine" , "aquamarine", "darkslateblue", "darkgreen" , "magenta", "purple" , "green", "brown") 

#############################


phes= c("Lutein2",
        "Lutein1",
        "ACrypto",
        "BCrypto",
        "alphabeta2",
        "alphabeta1",
        "LutRB",
        "ZeaRB",
        "ZeaBC",
        "TotBC",
        "GAE",
        "TKW22",
        "TKW23")

library(readxl)

#excel file of QTL information
BioQTL <- read_excel("~/Downloads/BioQTL.xlsx", 
                     sheet = "Sheet2")


#plot basic QTL map for ex) TKW 
BioQTLtkw<- subset(BioQTL, phenotype =="TKW"|
                     phenotype== "TKW23")

with(BioQTLtkw, segmentsOnMap(qtlfile, phe = phenotype, chr = chr,
                              l = lowposition, h = highposition, legendCex = 0.8,
                              peakcM = pos,
                              tick.width = .1,  chrBuffer = c(.15,.2), col = "red", legendPosition = "bottomright", lwd= 4))


View(BioQTL)
# segmentsOnMap(cross=qtlfile$geno, phe=BioQTL$phenotype,
#               col = col, ##add more colors if more phenotypes
#               chr=BioQTL$chr, 
#               l = BioQTL$lowposition, 
#               lwd=5,
#               chrBuffer = c(0.2, 0.2),
#               h = BioQTL$highposition,
#               legendPosition = "bottomright",
#               legendCex= 1,
#               leg.lwd = 5,
#               leg.inset=.05)

############################################ PCA of SAP data, PCAs from GAPIT
library(readr)
GAPIT_Genotype_PCA <- read_csv("~/Desktop/GWAS/GAPIT.Genotype.PCA.csv")
View(GAPIT_Genotype_PCA)

colnames(GAPIT_Genotype_PCA)[1]= "PI"

#PCAs from GAPIT data

pheno <- left_join(GAPIT_Genotype_PCA, KSURelativeBio_average, by = "PI")
pheno$pheno<- ifelse(pheno$Lutein > 0, "yes", "no")
ggplot(pheno,
       aes(
         PC1, 
         PC2,
         color= pheno
       ))+
  geom_point(alpha=0.8) +
  scale_color_manual(values = c("yes" = "orange", "treated" = "grey")) + theme_minimal()




pheno2$Color <-as.factor(pheno2$Color)
ggplot(pheno,
       aes(PC1,PC3,
           color=Color,))+
  geom_point(alpha=0.8, size=3) + scale_colour_manual(values = peri_clr_vec, breaks=c("brown","red","white", "yellow")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

################################## LOD score heatmap
####### uses scanone output from R/QTL

s2 <- tibble::rownames_to_column(s1, "BPpos")
s1$Bpos <- sub(".*_", "", s2$BPpos)

colnames(s1)[1]="Chr"
s1$Bpos= as.numeric(s1$Bpos)

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(scales)

plot_lod_heatmap <- function(s1, x_var = c("pos", "Bpos")) {
  x_var <- match.arg(x_var)
  x_sym <- rlang::sym(x_var)
  
  stopifnot(all(c("Chr", "pos", "Bpos") %in% names(s1)))
  phenotype_cols <- setdiff(names(s1), c("Chr", "pos", "Bpos"))
  if (!length(phenotype_cols)) stop("No phenotype columns beyond Chr/pos/Bpos found.")
  
  # Long format
  long_df <- s1 %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(dplyr::all_of(phenotype_cols),
                        names_to = "Phenotype", values_to = "LOD") %>%
    dplyr::mutate(
      Chr = factor(Chr, levels = unique(Chr)),
      Phenotype = factor(Phenotype, levels = phenotype_cols),
      x = if (x_var == "Bpos") !!x_sym / 1e6 else !!x_sym   # bp → Mb
    )
  
  # Tile widths
  tile_widths <- long_df %>%
    dplyr::distinct(Chr, x) %>%
    dplyr::arrange(Chr, x) %>%
    dplyr::group_by(Chr) %>%
    dplyr::summarise(
      dx = {
        diffs <- diff(x)
        if (length(diffs) >= 1 && any(is.finite(diffs))) {
          md <- stats::median(diffs[is.finite(diffs)], na.rm = TRUE)
          if (is.finite(md) && md > 0) md else (max(x, na.rm=TRUE)-min(x, na.rm=TRUE))*0.01
        } else {
          (max(x, na.rm=TRUE)-min(x, na.rm=TRUE))*0.01
        }
      },
      xmax = max(x, na.rm = TRUE),
      .groups = "drop"
    )
  
  long_df <- dplyr::left_join(long_df, tile_widths, by = "Chr")
  
  # Plot
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = x, y = Phenotype, fill = LOD)) +
    ggplot2::geom_tile(height = 0.9, width = pmax(long_df$dx, .Machine$double.eps)) +
    ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient(
      low = "yellow",
      high = "blue",
      na.value = "grey90"
    ) +
    ggplot2::labs(
      x = if (x_var == "pos") "Position (cM)" else "Position (Mb)",
      y = "Phenotype",
      fill = "LOD",
      title = "Genome-wide LOD heatmap by phenotype (rows) and chromosome (columns)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  # Format Mb axis nicely
  if (x_var == "Bpos") {
    p <- p + ggplot2::scale_x_continuous(
      labels = scales::label_number(accuracy = 1, suffix = "")
    )
  }
  
  # Add vertical black line at the right edge of each chromosome
  p <- p + ggplot2::geom_vline(data = tile_widths, 
                               ggplot2::aes(xintercept = xmax), 
                               color = "black", linewidth = 0.3, inherit.aes = FALSE)
  
  return(p)
}


plot_lod_heatmap(s1, x_var = "Bpos")

########################## LOD score heat map with labeled QTL and monomorphic markers at the bottom
```{r monomorphic markers}
diffm<- s1[,c(1,34)]
# 
Allmarkers <- read.csv("~/Desktop/qtl/Allmarkers.csv")
# View(Allmarkers)
diffm$Chr= as.numeric(diffm$Chr)
# 
colnames(Allmarkers)[3]= "Bpos"
unique_df1 <- anti_join(Allmarkers , diffm, by = c("Chr", "Bpos"))
BioKaryo <- read.csv("~/Downloads/BioKaryo.csv")
#View(BioKaryo)
BIOQTL <- read.csv("~/Desktop/qtl/BIOQTL.csv")

# Dependencies
library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(scales)
library(patchwork)  
#install.packages("patchwork") 

# s1: wide table with columns Chr, pos (cM), Bpos (bp), and phenotype columns (LOD values)
# BIOQTL: data.frame with columns phenotype, LOD, est, chr, Peak (bp)
# BioKaryo: data.frame with columns Chr, Start, End (bp)  # CE_* ignored
# unique_df1: data.frame with columns Chr, Bpos (bp) -> markers shown as grey ticks
plot_lod_heatmap <- function(s1,
                             x_var = c("pos", "Bpos"),
                             BIOQTL = NULL,
                             BioKaryo = NULL,
                             unique_df1 = NULL) {
  x_var <- match.arg(x_var)
  x_sym <- rlang::sym(x_var)
  
  stopifnot(all(c("Chr", "pos", "Bpos") %in% names(s1)))
  phenotype_cols <- setdiff(names(s1), c("Chr", "pos", "Bpos"))
  if (!length(phenotype_cols)) stop("No phenotype columns beyond Chr/pos/Bpos found.")
  
  # ---- Long format for heatmap ----
  long_df <- s1 %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(phenotype_cols),
      names_to = "Phenotype",
      values_to = "LOD"
    ) %>%
    dplyr::mutate(
      Chr = factor(Chr, levels = unique(Chr)),
      Phenotype = factor(Phenotype, levels = phenotype_cols),
      x = if (x_var == "Bpos") !!x_sym / 1e6 else !!x_sym  # bp → Mb when needed
    )
  
  # Tile width per chromosome (median spacing; fallback 1% span)
  tile_widths <- long_df %>%
    dplyr::distinct(Chr, x) %>%
    dplyr::arrange(Chr, x) %>%
    dplyr::group_by(Chr) %>%
    dplyr::summarise(
      dx = {
        diffs <- diff(x)
        if (length(diffs) >= 1 && any(is.finite(diffs))) {
          md <- stats::median(diffs[is.finite(diffs)], na.rm = TRUE)
          if (is.finite(md) && md > 0) md else (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) * 0.01
        } else {
          (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) * 0.01
        }
      },
      .groups = "drop"
    )
  long_df <- dplyr::left_join(long_df, tile_widths, by = "Chr")
  
  # ---- Karyotype boundaries (BioKaryo) ----
  karyo <- NULL
  if (!is.null(BioKaryo) && x_var == "Bpos") {
    required_k <- c("Chr", "Start", "End")
    if (!all(required_k %in% names(BioKaryo))) {
      stop("BioKaryo must contain columns: Chr, Start, End (bp).")
    }
    obs_ranges <- long_df %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(minx = min(x, na.rm = TRUE), maxx = max(x, na.rm = TRUE), .groups = "drop")
    
    karyo <- BioKaryo %>%
      dplyr::mutate(
        Chr = factor(Chr, levels = levels(long_df$Chr)),
        x_start = Start / 1e6, x_end = End / 1e6
      ) %>%
      dplyr::filter(!is.na(Chr)) %>%
      dplyr::left_join(obs_ranges, by = "Chr") %>%
      dplyr::mutate(
        x_start = dplyr::coalesce(x_start, minx),
        x_end   = dplyr::coalesce(x_end,   maxx)
      )
  }
  
  # --------- TOP: Heatmap (keep labeled x-axis here) ----------
  p_top <- ggplot2::ggplot(long_df, ggplot2::aes(x = x, y = Phenotype, fill = LOD)) +
    ggplot2::geom_tile(height = 0.9, width = pmax(long_df$dx, .Machine$double.eps)) +
    ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient(low = "yellow", high = "blue", na.value = "grey90") +
    ggplot2::labs(
      x = if (x_var == "pos") "Position (cM)" else "Position (Mb)",
      y = "Phenotype",
      fill = "LOD",
      title = "Genome-wide LOD heatmap with QTL markers"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(b = 4)  # small gap above the marker strip
    )
  if (x_var == "Bpos") {
    p_top <- p_top + ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 1, suffix = ""))
  }
  if (!is.null(karyo)) {
    p_top <- p_top + ggplot2::geom_rect(
      data = karyo,
      ggplot2::aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf),
      fill = NA, color = NA, inherit.aes = FALSE
    )
  }
  # Vertical separator at right edge of each chromosome
  vlines_df <- if (!is.null(karyo)) {
    karyo %>% dplyr::transmute(Chr, xmax = x_end)
  } else {
    long_df %>% dplyr::group_by(Chr) %>%
      dplyr::summarise(xmax = max(x, na.rm = TRUE), .groups = "drop")
  }
  p_top <- p_top + ggplot2::geom_vline(
    data = vlines_df, ggplot2::aes(xintercept = xmax),
    color = "black", linewidth = 0.3
  )
  
  # QTL overlay (option B: outline triangles, color only)
  if (!is.null(BIOQTL)) {
    required_q <- c("phenotype", "LOD", "est", "chr", "Peak")
    if (!all(required_q %in% names(BIOQTL))) {
      stop("BIOQTL must contain columns: phenotype, LOD, est, chr, Peak")
    }
    qtl_df <- BIOQTL %>%
      dplyr::mutate(
        Chr = factor(chr, levels = levels(long_df$Chr)),
        Phenotype = factor(phenotype, levels = levels(long_df$Phenotype)),
        xqtl = if (x_var == "Bpos") Peak / 1e6 else NA_real_,
        direction = dplyr::case_when(est > 0 ~ "pos", est < 0 ~ "neg", TRUE ~ NA_character_)
      ) %>%
      dplyr::filter(!is.na(Chr), !is.na(Phenotype), !is.na(direction), !is.na(xqtl))
    if (x_var == "pos") {
      warning("BIOQTL$Peak is in bp; supply cM peak positions to overlay when x_var = 'pos'. Skipping QTL overlay.")
    } else {
      p_top <- p_top +
        ggplot2::geom_point(
          data = qtl_df,
          ggplot2::aes(x = xqtl, y = Phenotype, shape = direction, color = direction),
          size = 2.6, stroke = 1.0, inherit.aes = FALSE
        ) +
        ggplot2::scale_shape_manual(values = c(pos = 2, neg = 6), guide = "none") +
        ggplot2::scale_color_manual(values = c(pos = "green3", neg = "red3"), guide = "none")
    }
  }
  
  # --------- BOTTOM: Marker strip as grey tick marks ----------
  p_bottom <- NULL
  if (!is.null(unique_df1)) {
    markers <- unique_df1 %>%
      dplyr::transmute(
        Chr = factor(Chr, levels = levels(long_df$Chr)),
        x_m = if (x_var == "Bpos") Bpos / 1e6 else NA_real_
      ) %>%
      dplyr::filter(!is.na(Chr), is.finite(x_m))
    
    if (nrow(markers) > 0) {
      p_bottom <- ggplot2::ggplot(markers, ggplot2::aes(x = x_m, y = 0)) +
        # short vertical tick marks at each marker position
        ggplot2::geom_segment(aes(xend = x_m, yend = 0.6), color = "grey40", linewidth = 0.4) +
        ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
        ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = margin(t = 0)  # sits tight beneath the top plot's x-axis labels
        )
      
      if (!is.null(karyo)) {
        p_bottom <- p_bottom + ggplot2::geom_rect(
          data = karyo,
          ggplot2::aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf),
          fill = NA, color = NA, inherit.aes = FALSE
        )
      }
      # (No x labels here; labels are on the top plot per your request)
      
      # Optional: same end-of-chromosome separator (kept off to avoid clutter)
      # p_bottom <- p_bottom + ggplot2::geom_vline(
      #   data = vlines_df, ggplot2::aes(xintercept = xmax),
      #   color = "black", linewidth = 0.3
      # )
    } else if (x_var == "pos") {
      warning("unique_df1 only has Bpos (bp); supply cM marker positions to show markers when x_var = 'pos'.")
    }
  }
  
  # Combine: top plot retains labeled x-axis; marker strip sits directly beneath it
  if (!is.null(p_bottom)) {
    return(p_top / p_bottom + patchwork::plot_layout(heights = c(5, 0.2)))
  } else {
    return(p_top)
  }
}

# Example:
p <- plot_lod_heatmap(
  s1,
  x_var = "Bpos",
  BIOQTL = BIOQTL,
  BioKaryo = BioKaryo,
  unique_df1 = unique_df1
)
print(p)

png("LODheatmapv4.png", width = 4000, height = 2000, res = 300) # width/height in pixels, res in DPI

plot(p, main = "hQTL")

# Close the device to save the plot
dev.off()




