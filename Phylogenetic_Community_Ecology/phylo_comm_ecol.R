## loading libraries ===========
if(!require("xfun")) install.packages("xfun")
options(repos = c(
  rtrees = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
xfun::pkg_attach2(c("tidyverse", "picante", "ape", 
                    "rtrees", # to get phylogeny
                    "hillR", # calculate hill number based diversity values
                    "phyr"
))

# loading data =========

veg_long = read.csv("data/veg_long.csv")
veg_wide = read.csv("data/veg_wide.csv", row.names = 1)
traits = read.csv("data/traits.csv")
envi = read.csv("data/envi.csv")


head(veg_long)
head(veg_wide)
dim(veg_wide) # 30 sites, 55 species 
head(traits)
head(envi)

# things to consider:
# - all species have traits
# - all sites have environmental data
# - vegetation data do not include species or sites with no observation

# here are the species
sort(colnames(veg_wide))

# get a phylogeny for these species
# note that there are species like: Carex_spp
phy = rtrees::get_tree(sp_list = sort(colnames(veg_wide)), taxon = "plant")
# which ones are grafted?
phy$graft_status
plot(phy, cex = 0.6)
plot(phy, type = "fan", cex = 0.6)
# all species should be in the phylogeny
all(colnames(veg_wide) %in% phy$tip.label)
# to calculate phylogenetic diversity or phylogenetic signals, phylogenies derived
# in this way works well: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2788

# Calculate diversity for each community =====
## species diversity
hill_taxa(comm = veg_wide, q = 0) # richness
hill_taxa(comm = veg_wide, q = 1) # shannon
hill_taxa(comm = veg_wide, q = 2) # simpson
# put them in a data frame
div = data.frame(
  site = rownames(veg_wide),
  sp_richness = hill_taxa(comm = veg_wide, q = 0),
  sp_shannon = hill_taxa(comm = veg_wide, q = 1),
  sp_simpson = hill_taxa(comm = veg_wide, q = 2)
)

## functional diversity
# hill_func(comm = veg_wide, traits = traits)
# error above?!
# we need to put species names as row names ...
row.names(traits) = traits$sp
str(traits)

hf0 = hill_func(comm = veg_wide, traits = traits[, -1], q = 0)
hf0
hf1 = hill_func(comm = veg_wide, traits = traits[, -1], q = 1)
hf2 = hill_func(comm = veg_wide, traits = traits[, -1], q = 2)

div$FDis = hf0["FDis",]
div$FD_0 = hf0["FD_q",]
div$FD_1 = hf1["FD_q",]
div$FD_2 = hf2["FD_q",]

## phylogenetic diversity
hill_phylo(comm = veg_wide, tree = phy, q = 0)

div$PD_0 = hill_phylo(comm = veg_wide, tree = phy, q = 0)
div$PD_1 = hill_phylo(comm = veg_wide, tree = phy, q = 1)
div$PD_2 = hill_phylo(comm = veg_wide, tree = phy, q = 2)

# other common phylogenetic diversity metrics
picante::psv(samp = veg_wide, tree = phy)
div$psv = picante::psv(samp = veg_wide, tree = phy)$PSVs
div$pse = picante::pse(samp = veg_wide, tree = phy)$PSEs
div$mpd_abund = picante::mpd(samp = veg_wide, dis = cophenetic(phy), abundance.weighted = TRUE)
div$mpd_pa = picante::mpd(samp = veg_wide, dis = cophenetic(phy), abundance.weighted = FALSE)
div$mntd_abund = picante::mntd(samp = veg_wide, dis = cophenetic(phy), abundance.weighted = TRUE)
div$mntd_pa = picante::mntd(samp = veg_wide, dis = cophenetic(phy), abundance.weighted = FALSE)

# explore data
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  z = na.omit(data.frame(x = x, y = y))
  r <- abs(cor(z$x, z$y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = 4 * cex.cor * r)
}

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
}

pairs(div[, -1], lower.panel = panel.cor, diag.panel = panel.hist)

# join diversity values with environmental variables ====
div_envi = left_join(div, envi, by = "site")
div_envi
div_envi = as_tibble(div_envi)
div_envi

# one can then explore the relationships between diversity and environmental variables 
# not going to do here ...


# PGLMMs ====
veg_long = mutate(veg_long, freq_log = log(freq + 1))

if(!file.exists("data/m1.rds")){
  m1 <- pglmm(freq_log ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), data = veg_long, cov_ranef = list(sp = phy))
  saveRDS(m1, "data/m1.rds")
} else {
  m1 = readRDS("data/ms.rds")
}

summary(m1)
# to test the significance of attraction
pglmm_profile_LRT(m1, re.number = 4)

# so there is strong attraction in these communities
veg_long = left_join(veg_long, traits, by = "sp")

if(!file.exists("data/m1.rds")){
  m2 <- pglmm(freq_log ~ 1 + circ + l.width + l.thick + animal.disp +
                (1|sp__) + (1|site) + (1|sp__@site) +
                (circ|site) + (l.width|site) + (l.thick|site) + (animal.disp|site), 
              data = veg_long, cov_ranef = list(sp = phy))
  saveRDS(m2, "data/m2.rds")
} else {
  m2 = readRDS("data/m2.rds")
}

summary(m2)
# to test the significance of attraction
pglmm_profile_LRT(m2, re.number = 4)

# to learn more about pglmm(), see the documentation with ?pglmm

