# R Script: Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: phylogenetic comparative evidence from scorpion pincers
# Stenio I A Foerster
# https://foersterst.github.io
# May 2026


# Libraries & Data ----------------------------------------------------------------------------

library(marginaleffects)
library(bayestestR)
library(tidyverse)
library(phytools)
library(flextable)
library(patchwork)
library(rptR)
library(ggdist)
library(ggsci)
library(coda)
library(brms)

# read in the tree
tree <- read.tree("data/scorp-complete-tree.tre") # with imputation
# tree <- read.tree("data/scorp-complete-noimp.tre") # without imputation

# read in the trait data
traits_data <- read_csv("data/scorp-complete.csv") # with imputation
# traits_data <- read_csv("data/scorp-complete-noimp.csv") # without imputation
traits_data <- traits_data[match(tree$tip.label, traits_data$species), ]
identical(traits_data$species, tree$tip.label)

# set theme for all graphs
set_theme(theme_classic(base_line_size = 0.3, base_rect_size = 0.3))
theme_update(
  axis.line = element_blank(),
  panel.background = element_rect(fill = "#eae3d2"),
  panel.grid.major = element_line(color = "white", linewidth = 0.2),
  axis.title = element_text(color = "black", size = 9),
  axis.text = element_text(color = "black", size = 8.5),
  legend.title = element_text(color = "black", size = 9),
  legend.text = element_text(color = "black", size = 8.5),
  legend.key.height = unit(3.6, "mm"),
  legend.key.width = unit(5, "mm"),
  legend.key.spacing.y = unit(1, "mm")
)


# 1. RevBayes Input ---------------------------------------------------------------------------

# Make input files for running OU models in RevBayes

# calculate phylogenetic residuals
# traits_data$cheL_resid <- phyl.resid(tree = tree, x = setNames(traits_data$log10_carL, traits_data$species), Y = setNames(traits_data$log10_cheL, traits_data$species), method = "lambda")$resid[, 1]
# traits_data$cheW_resid <- phyl.resid(tree = tree, x = setNames(traits_data$log10_carL, traits_data$species), Y = setNames(traits_data$log10_cheW, traits_data$species), method = "lambda")$resid[, 1]
# traits_data$cheD_resid <- phyl.resid(tree = tree, x = setNames(traits_data$log10_carL, traits_data$species), Y = setNames(traits_data$log10_cheD, traits_data$species), method = "lambda")$resid[, 1]

# RevBayes requires trait data in nexus format. I will create a backbone nexus file with species and traits, and then add the nexus notation manually with a text editor.

# traits_data |>
#   select(species, cheL_resid, cheW_resid, cheD_resid) |>
#   #write.table(file = "rb/data/scorp-chela.nex", quote = F, sep = "\t", row.names = F, col.names = F) # with imputation
#   #write.table(file = "rb/data/scorp-chela-noimp.nex", quote = F, sep = "\t", row.names = F, col.names = F) # without imputation

# export tree in nexus format
# write.nexus(tree, file = "rb/data/scorp-tree.nex") # with imputation
# write.nexus(tree, file = "rb/data/scorp-tree-noimp.nex") # without imputation


# 2. RevBayes Post-process --------------------------------------------------------------------

# summarize RevBayes results using graphs and tables

# file list
tt <- c("cheL", "cheW", "cheD")

# tab to store results
mcmc_tab <- data.frame()

for (i in 1:length(tt)) {
  # list & read rb files
  f1 <- list.files("rb/output/", pattern = tt[i], full.names = T)
  f1 <- f1[!grepl("noimp", f1)] # with imputation
  # f1 <- f1[grepl("noimp", f1)] # without imputation

  f2 <- f1[grepl("run", f1)]
  f1 <- f1[!grepl("run", f1)]
  x1 <- read.table(f1, header = T)[, -c(1, 2, 12)]
  x2 <- read.table(f2[1], header = T)[, -c(1, 11)]
  x3 <- read.table(f2[2], header = T)[, -c(1, 11)]

  # temporary data frame to store results
  tmp <- as.data.frame(HPDinterval(mcmc(x1)))
  tmp$median <- apply(x1, MARGIN = 2, median, na.rm = T) # median
  tmp$ess <- effectiveSize(mcmc(x1)) # effective sample size
  tmp <- data.frame(trait = tt[i], param = rownames(tmp), tmp, gelman.diag(mcmc.list(mcmc(x2), mcmc(x3)))$psrf) # Gelman & Rubin's convergence diagnostic
  rownames(tmp) <- NULL
  # add the results
  mcmc_tab <- rbind(mcmc_tab, tmp)
}

rm(tt, f1, f2, x1, x2, x3, tmp, i)
mcmc_tab

# export for supp
# write_csv(mcmc_tab, file = "supp/mcmc-tab.csv")
# write_csv(mcmc_tab, file = "supp/mcmc-tab-noimp.csv")


## 2.1. Tests ----

# read in the combined mcmc runs
cheL_mcmc <- read_table("rb/output/cheL.log") # with imputation
cheW_mcmc <- read_table("rb/output/cheW.log") # with imputation
cheD_mcmc <- read_table("rb/output/cheD.log") # with imputation

# cheL_mcmc <- read_table("rb/output/cheL-noimp.log") # without imputation
# cheW_mcmc <- read_table("rb/output/cheW-noimp.log") # without imputation
# cheD_mcmc <- read_table("rb/output/cheD-noimp.log") # without imputation

# NB: P-values are frequentist-analogue, two-sided P values

# sigma2
p_direction(c(cheD_mcmc$sigma - cheW_mcmc$sigma), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$sigma - cheL_mcmc$sigma), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$sigma - cheL_mcmc$sigma), as_p = T, null = 0) # width vs length

# alpha
p_direction(c(cheD_mcmc$alpha - cheW_mcmc$alpha), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$alpha - cheL_mcmc$alpha), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$alpha - cheL_mcmc$alpha), as_p = T, null = 0) # width vs length

# rho
p_direction(c(cheD_mcmc$rho - cheW_mcmc$rho), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$rho - cheL_mcmc$rho), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$rho - cheL_mcmc$rho), as_p = T, null = 0) # width vs length

# phylogenetic half-life
p_direction(c(cheD_mcmc$t_half - cheW_mcmc$t_half), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$t_half - cheL_mcmc$t_half), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$t_half - cheL_mcmc$t_half), as_p = T, null = 0) # width vs length

# eta
p_direction(c(cheD_mcmc$eta - cheW_mcmc$eta), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$eta - cheL_mcmc$eta), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$eta - cheL_mcmc$eta), as_p = T, null = 0) # width vs length

# expected variance (V)
p_direction(c(cheD_mcmc$V - cheW_mcmc$V), as_p = T, null = 0) # height vs width
p_direction(c(cheD_mcmc$V - cheL_mcmc$V), as_p = T, null = 0) # height vs length
p_direction(c(cheW_mcmc$V - cheL_mcmc$V), as_p = T, null = 0) # width vs length

# median absolute contrasts
median(c(cheD_mcmc$eta - cheW_mcmc$eta))
median(c(cheD_mcmc$V - cheW_mcmc$V))
median(c(cheD_mcmc$alpha - cheW_mcmc$alpha))
median(c(cheD_mcmc$sigma - cheW_mcmc$sigma))

## 2.2. Graphs ----

# read in combined (post-burnin) runs
x1 <- read_table("rb/output/cheL.log") # with imputation
x2 <- read_table("rb/output/cheW.log") # with imputation
x3 <- read_table("rb/output/cheD.log") # with imputation

# x1 <- read_table("rb/output/cheL-noimp.log") # without imputation
# x2 <- read_table("rb/output/cheW-noimp.log") # without imputation
# x3 <- read_table("rb/output/cheD-noimp.log") # without imputation

# add trait tag
x1$Trait <- "cheL"
x2$Trait <- "cheW"
x3$Trait <- "cheD"

# combine & clean
mcmc_runs <- rbind(x1, x2, x3)
rm(x1, x2, x3)

# some editing
mcmc_runs |>
  mutate(
    Trait = factor(
      recode(
        Trait,
        "cheL" = "Length",
        "cheW" = "Width",
        "cheD" = "Height"
      ),
      levels = c("Length", "Width", "Height")
    )
  ) -> mcmc_runs

# sigma2
ggplot(mcmc_runs, aes(x = sigma)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  scale_color_jama() +
  scale_fill_jama() +
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$sigma), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$sigma), y = -0.03, label = "b", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$sigma), y = -0.38, label = "c", size = 3.2) -> plot_sigma

# stationary variance (eta)
ggplot(mcmc_runs, aes(x = eta)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  scale_color_jama() +
  scale_fill_jama() +
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$eta), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$eta), y = -0.03, label = "a", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$eta), y = -0.38, label = "b", size = 3.2) -> plot_eta

# expected tip variance (V)
ggplot(mcmc_runs, aes(x = V)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  scale_color_jama() +
  scale_fill_jama() +
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$V), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$V), y = -0.03, label = "a", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$V), y = -0.38, label = "b", size = 3.2) -> plot_V

# alpha
ggplot(mcmc_runs, aes(x = alpha)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  # with imputation
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$alpha), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$alpha), y = -0.03, label = "a", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$alpha), y = -0.38, label = "b", size = 3.2) +
  # without imputation
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$alpha), y = 0.3, label = "a", size = 3.2) + # height
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$alpha), y = -0.03, label = "b", size = 3.2) + # width
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$alpha), y = -0.38, label = "a", size = 3.2) +
  scale_color_jama() +
  scale_fill_jama() -> plot_alpha

# rho
ggplot(mcmc_runs, aes(x = rho)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  # with imputation
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$rho), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$rho), y = -0.03, label = "a", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$rho), y = -0.38, label = "b", size = 3.2) +
  # without imputation
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$rho), y = 0.3, label = "a", size = 3.2) + # height
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$rho), y = -0.03, label = "b", size = 3.2) + # width
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$rho), y = -0.38, label = "a", size = 3.2) +
  scale_color_jama() +
  scale_fill_jama() -> plot_rho

# half-life
ggplot(mcmc_runs, aes(x = t_half)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  # with imputation
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$t_half), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$t_half), y = -0.03, label = "a", size = 3.2) + # width
  annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$t_half), y = -0.38, label = "b", size = 3.2) +
  # without imputation
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Height", ]$t_half), y = 0.3, label = "a", size = 3.2) + # height
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Width", ]$t_half), y = -0.03, label = "b", size = 3.2) + # width
  # annotate(geom = "text", x = median(mcmc_runs[mcmc_runs$Trait == "Length", ]$t_half), y = -0.38, label = "a", size = 3.2) +
  scale_color_jama() +
  scale_fill_jama() -> plot_half

# OU panel
plot_sigma <- plot_sigma + labs(y = NULL, x = expression("Evolutionary rate" ~ (sigma^2)))
plot_alpha <- plot_alpha + labs(y = NULL, x = expression("Stabilizing selection" ~ (alpha)))
plot_eta <- plot_eta + labs(y = NULL, x = expression("Stationary variance" ~ (eta)))
plot_V <- plot_V + labs(y = NULL, x = "Expected variance (V)")
plot_half <- plot_half + labs(y = NULL, x = expression("Phylogenetic half-life" ~ (t[1 / 2])))
plot_rho <- plot_rho + labs(y = NULL, x = expression("Percent decr. in trait var." ~ (rho)))

(plot_sigma + plot_alpha + plot_eta + plot_V + plot_half + plot_rho) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) -> pp1

# export with imputation
# ggsave(filename = "output/ou-params.png", plot = pp1, width = 15, height = 17, units = "cm", dpi = 900)
# ggsave(filename = "output/ou-params.pdf", plot = pp1, width = 15, height = 17, units = "cm", dpi = 900)
# export without imputation
# ggsave(filename = "output/ou-params-noimp.png", plot = pp1, width = 15, height = 17, units = "cm", dpi = 900)
# ggsave(filename = "output/ou-params-noimp.pdf", plot = pp1, width = 15, height = 17, units = "cm", dpi = 900)


# 3. Evolutionary Allometry -------------------------------------------------------------------

# Bayesian multivaraite allometric model

# see info in the brms vignette: https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html

# covariance matrix of species
A <- ape::vcv(tree, corr = TRUE)

# copy of the data
df <- traits_data

# center carL to improve mixing and to ensure that phylogenetic random effects are evaluated at mean size (not at zero)
# note that scale = FALSE, so nothing changes about the biological meaning of the slope (unit)
df$log10_carL <- scale(df$log10_carL, center = TRUE, scale = FALSE)

# fit the model
brm_fit <- brm(
  mvbind(log10_cheL, log10_cheW, log10_cheD) ~ log10_carL +
    (1 | gr(species, cov = A)),
  data = df,
  data2 = list(A = A),
  family = gaussian(),
  seed = 123,
  chains = 2,
  warmup = 2000,
  iter = 4000,
  cores = 4
)

## 3.1. Diagnostics ----

# Rhat and effective sample size (ESS)
summary(brm_fit)

# trace plots
bayesplot::mcmc_trace(
  ps,
  pars = colnames(ps)[grepl("b_|rescor", colnames(ps))],
  facet_args = list(ncol = 3, labeller = label_parsed)
) +
  facet_text(size = 7) -> tp

# save
# ggsave(filename = "supp/brms-trace.png", plot = tp, width = 16, height = 15, units = "cm", dpi = 900)
# ggsave(filename = "supp/brms-trace.pdf", plot = tp, width = 16, height = 15, units = "cm", dpi = 900)

# posterior predictive checks
pr1 <- pp_check(brm_fit, resp = "log10cheL", ndraws = 100)
pr2 <- pp_check(brm_fit, resp = "log10cheW", ndraws = 100)
pr3 <- pp_check(brm_fit, resp = "log10cheD", ndraws = 100)

pr1 <- pr1 + labs(x = expression(log[10] ~ "chela length"))
pr2 <- pr2 + labs(x = expression(log[10] ~ "chela width"))
pr3 <- pr3 + labs(x = expression(log[10] ~ "chela height"))

(pr1 + pr2 + pr3) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) -> bpp1

# save
# ggsave(filename = "supp/brms-slopes-ppc.png", plot = bpp1, width = 16, height = 7, units = "cm", dpi = 900)
# ggsave(filename = "supp/brms-slopes-ppc.pdf", plot = bpp1, width = 16, height = 7, units = "cm", dpi = 900)

# export summary table to be shown in the supp
# j <- summary(brm_fit)$fixed
# j$param <- rownames(j)
# write_csv(j, file = "supp/brms-slopes.csv") # with imputation
# write_csv(j, file = "supp/brms-slopes-noimp.csv") # without imputation


## 3.2. Tests ----

# get the posterior distributions
ps <- as_tibble(as_draws_df(brm_fit))

# NB: P-values are frequentist-analogue, two-sided P values

# slopes
p_direction(c(ps$b_log10cheD_log10_carL - ps$b_log10cheW_log10_carL), as_p = T, null = 0) # height vs width
p_direction(c(ps$b_log10cheD_log10_carL - ps$b_log10cheL_log10_carL), as_p = T, null = 0) # height vs length
p_direction(c(ps$b_log10cheW_log10_carL - ps$b_log10cheL_log10_carL), as_p = T, null = 0) # width vs length

# absolute contrasts
summary(c(ps$b_log10cheD_log10_carL - ps$b_log10cheW_log10_carL))
summary(c(ps$b_log10cheD_log10_carL - ps$b_log10cheL_log10_carL))
summary(c(ps$b_log10cheW_log10_carL - ps$b_log10cheL_log10_carL))

# summary
posterior_summary(ps$b_log10cheD_log10_carL)
posterior_summary(ps$b_log10cheW_log10_carL)
posterior_summary(ps$b_log10cheL_log10_carL)

# residual correlations
p_direction(c(ps$rescor__log10cheL__log10cheD - ps$rescor__log10cheL__log10cheW), as_p = T, null = 0) # LH vs LW
p_direction(c(ps$rescor__log10cheL__log10cheD - ps$rescor__log10cheW__log10cheD), as_p = T, null = 0) # LH vs WH
p_direction(c(ps$rescor__log10cheL__log10cheW - ps$rescor__log10cheW__log10cheD), as_p = T, null = 0) # LW vs WH

# summary
posterior_summary(ps$rescor__log10cheL__log10cheD)
posterior_summary(ps$rescor__log10cheL__log10cheW)
posterior_summary(ps$rescor__log10cheW__log10cheD)


## 3.3. Graphs ----

# posterior slopes #

ps |>
  dplyr::select(b_log10cheL_log10_carL, b_log10cheW_log10_carL, b_log10cheD_log10_carL) |>
  rename(
    "Length" = b_log10cheL_log10_carL,
    "Width" = b_log10cheW_log10_carL,
    "Height" = b_log10cheD_log10_carL
  ) |>
  mutate(type = "Slope") |>
  pivot_longer(!type, names_to = "Trait", values_to = "posterior") |>
  mutate(Trait = factor(Trait, levels = c("Length", "Width", "Height"))) |>
  ggplot(aes(x = posterior)) +
  stat_halfeye(aes(fill = Trait, color = Trait),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  # with imputation
  annotate(geom = "text", x = median(ps$b_log10cheD_log10_carL), y = 0.3, label = "a", size = 3.2) + # height
  annotate(geom = "text", x = median(ps$b_log10cheL_log10_carL), y = -0.38, label = "b", size = 3.2) + # length
  annotate(geom = "text", x = median(ps$b_log10cheW_log10_carL), y = -0.03, label = "b", size = 3.2) +
  # without imputation
  # annotate(geom = "text", x = median(ps$b_log10cheD_log10_carL), y = 0.3, label = "a", size = 3.2) + # height
  # annotate(geom = "text", x = median(ps$b_log10cheL_log10_carL), y = -0.38, label = "a", size = 3.2) + # length
  # annotate(geom = "text", x = median(ps$b_log10cheW_log10_carL), y = -0.03, label = "a", size = 3.2) +
  labs(y = NULL, x = "Allometric slope") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  scale_color_jama() +
  scale_fill_jama() -> pp2a

# residual correlation #

ps |>
  dplyr::select(
    rescor__log10cheL__log10cheD,
    rescor__log10cheL__log10cheW,
    rescor__log10cheW__log10cheD
  ) |>
  rename(
    "Length,Height" = rescor__log10cheL__log10cheD,
    "Length,Width" = rescor__log10cheL__log10cheW,
    "Width,Height" = rescor__log10cheW__log10cheD
  ) |>
  mutate(type = "rescorr") |>
  pivot_longer(!type, names_to = "Correlation", values_to = "posterior") |>
  ggplot(aes(x = posterior)) +
  stat_halfeye(aes(fill = Correlation, color = Correlation),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  # annotate(geom = "text", x = median(ps$rescor__log10cheL__log10cheD), y = -0.38, label = "b", size = 3.2) +
  # annotate(geom = "text", x = median(ps$rescor__log10cheL__log10cheW), y = -0.03, label = "b", size = 3.2) +
  # annotate(geom = "text", x = median(ps$rescor__log10cheW__log10cheD), y = 0.3, label = "a", size = 3.2) +
  labs(y = NULL, x = "Residual correlation") +
  scale_color_manual(values = c("#B24745FF", "#79AF97FF", "#6A6599FF")) +
  scale_fill_manual(values = c("#B24745FF", "#79AF97FF", "#6A6599FF")) -> pp2b

# combine graphs #

(pp2a + pp2b) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) -> pp2c

# save with imputation
# ggsave(filename = "output/allometric-slopes-multivar.png", plot = pp2c, width = 16, height = 7, units = "cm", dpi = 900)
# ggsave(filename = "output/allometric-slopes-multivar.pdf", plot = pp2c, width = 16, height = 7, units = "cm", dpi = 900)
# without imputation
# ggsave(filename = "output/allometric-slopes-multivar-noimp.png", plot = pp2c, width = 16, height = 7, units = "cm", dpi = 900)
# ggsave(filename = "output/allometric-slopes-multivar-noimp.pdf", plot = pp2c, width = 16, height = 7, units = "cm", dpi = 900)


# 4. Phylogenetic PCA -------------------------------------------------------------------------

# must be TRUE
identical(traits_data$species, tree$tip.label)

# input data
data_pca <- data.frame(
  cheL = traits_data$log10_cheL,
  cheW = traits_data$log10_cheW,
  cheD = traits_data$log10_cheD,
  carL = traits_data$log10_carL,
  row.names = traits_data$species
)

# ppca
pca_traits <- phyl.pca(tree = tree, Y = data_pca, method = "lambda")

# trait loadings
pca_traits$L

# export trait loadings to be presented in the supp
# rr <- data.frame(Trait = rownames(pca_traits$L), pca_traits$L)
# write_csv(rr, file = "supp/pca-loadings.csv")

# graph
factoextra::fviz_pca_biplot(
  X = as.princomp(pca_traits),
  geom = "point",
  arrow = 0.3,
  mean.point = F,
  col.ind = data_pca$carL,
  labelsize = 3,
  pointsize = 1.25,
  alpha.ind = 0.55,
  repel = T,
  col.var = "black",
  title = NULL
) +
  labs(x = "PC1 (90%)", y = "PC2 (6%)", color = expression(log[10] ~ carL)) +
  theme(
    axis.title = element_text(color = "black", size = 9),
    axis.text = element_text(color = "black", size = 8.5),
    legend.text = element_text(color = "black", size = 8.5),
    legend.title = element_text(color = "black", size = 9),
    panel.grid = element_blank()
    # panel.grid = element_line(linewidth = 0.25)
  ) +
  scale_color_viridis_c() -> pp3

# export
# ggsave(filename = "output/phylo-pca.png", plot = pp3, width = 13, height = 10, units = "cm", dpi = 900)
# ggsave(filename = "output/phylo-pca.pdf", plot = pp3, width = 13, height = 10, units = "cm", dpi = 900)


# 5. Repeatability --------------------------------------------------------


# Repeatability (R) tells whether species are good representations of true species values: the higher the better.
# An example from the literature is available here: https://doi.org/10.1111/evo.13865

# data males and females
dd <- read_csv("data/data-ssd.csv")

# total body size variance accounted for by differences between species after controlling for sex
rpt(
  carL ~ sex + (1 | species),
  grname = "species",
  data = dd,
  datatype = "Gaussian",
  nboot = 100
)

# total cheL variance accounted for by differences between species after controlling for body size and sex
rpt(
  cheL ~ carL + sex + (1 | species),
  grname = "species",
  data = dd,
  datatype = "Gaussian",
  nboot = 100
)

# total cheW variance accounted for by differences between species after controlling for body size and sex
rpt(
  cheW ~ carL + sex + (1 | species),
  grname = "species",
  data = dd,
  datatype = "Gaussian",
  nboot = 100
)

# total cheD variance accounted for by differences between species after controlling for body size and sex
rpt(
  cheD ~ carL + sex + (1 | species),
  grname = "species",
  data = dd,
  datatype = "Gaussian",
  nboot = 100
)


# END
