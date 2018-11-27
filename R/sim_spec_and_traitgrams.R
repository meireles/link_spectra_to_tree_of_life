library("phytools")
library("mvMORPH")
library("PEcAnRTM")

################################################################################
# Simulate Tree
################################################################################

# This seed outputs a reasonable tree shape
set.seed(114211)

# Base tree
tree_base = phytools::pbtree(n = 7, scale = 1)

picks     = c(1, 4)

tree_poly = phytools::paintSubTree(tree_base, node = picks[1], stem = TRUE,
                                   state = "orange", anc.state = "black")
tree_poly = phytools::paintSubTree(tree_poly, node = picks[2], stem = TRUE,
                                   state = "blue", anc.state = "black")
background_col  = rgb(0,0,0,0.4)
spec_colvec     = rep(background_col, Ntip(tree_poly))
spec_colvec[picks[1]]  = "orange"
spec_colvec[picks[2]]  = "blue"

# Unset seed
#set.seed(Sys.time())

################################################################################
# Simulate Traits
################################################################################

prospect_mu = c("N"   = 1.6,     # structural param. N of layers
                "Cab" = 20.0,    # chlorophyll a + b in ug cm^-2
                "Car" = 12.0,    # total carotenoids in ug cm^-2
                "Cw"  = 0.02,    # equivalent water thickness (EWT) in cm
                "Cm"  = 0.01)    # LMA in g cm^-2

sigma_f = c("N" = 0.3, "Cab" = 1, "Car" = 0.5, "Cw"  = 0.001, "Cm"  = 0.0005)
sigma_s = sigma_f
sigma_s[1] = sigma_s[1] / 2

# Sigma^2
sigma_f_sq = sigma_f * sigma_f
sigma_s_sq = sigma_s * sigma_s

# Alpha
alpha    = log(2)
alpha_s  = 1e-12
ou_alpha = matrix(data = diag( c(alpha, alpha_s, alpha_s, alpha_s, alpha_s)),
                  nrow = 5, ncol = 5,
                  dimnames = list(names(prospect_mu),
                                  names(prospect_mu)) )
## Params

params_bmf = list(ntraits      = 5,
                  sigma        = diag(sigma_f_sq),
                  theta        = prospect_mu,
                  names_traits = names(prospect_mu))

params_bms = list(ntraits      = 5,
                  sigma        = diag(sigma_s_sq),
                  theta        = prospect_mu,
                  names_traits = names(prospect_mu))

params_ouf = list(ntraits      = 5,
                  sigma        = diag(sigma_f_sq),
                  alpha        = ou_alpha,
                  theta        = prospect_mu,
                  names_traits = names(prospect_mu),
                  vcv          = "fixedRoot")

data_bmf = mvMORPH::mvSIM(tree  = tree_poly,
                          param = params_bmf,
                          model = "BM1",
                          nsim  = 1)

data_bms = mvMORPH::mvSIM(tree  = tree_poly,
                          param = params_bms,
                          model = "BM1",
                          nsim  = 1)

data_ouf = mvMORPH::mvSIM(tree  = tree_poly,
                          param = params_ouf,
                          model = "OU1",
                          nsim  = 1)

################################################################################
# Simulate Spectra
################################################################################

spec_bmf = apply(data_bmf, 1, function(x){
    PEcAnRTM::prospect(x, version = "5")[ , 1]
})

spec_bms = apply(data_bms, 1, function(x){
    PEcAnRTM::prospect(x, version = "5")[ , 1]
})

spec_ouf = apply(data_ouf, 1, function(x){
    PEcAnRTM::prospect(x, version = "5")[ , 1]
})

rownames(spec_bmf) = rownames(spec_bms) = rownames(spec_ouf) = 400:2500


########################################
# Pick Trait
#
# Only `N` has been simuated under an OUM model
# Others only have BM (fast & slow) and OU1
########################################

example_trait = setNames(object = "N", nm = "Leaf Structure (N)")

trait_plot    = example_trait
ylab          = names(example_trait)
ylim          = range(data_bmf[ , trait_plot],
                      data_bms[ , trait_plot],
                      data_ouf[ , trait_plot])
lwd           = 1.6
cex.main      = 1
fsize         = 0.02
spread        = FALSE

alpha_spec    = 0.2
ylim_spec     = range(range(spec_bmf),
                      range(spec_bms),
                      range(spec_ouf))

png("figures/sim_example_traitigram_and_spectra.png",
    width = 6 , height = 4, units = "in", res = 800)

par(mfcol = c(2, 3), mar = c(4, 5, 2, 2), oma = c(0, 0, 2, 0))

####################
phytools::phenogram(tree  = tree_poly,
                    x     = data_bmf[ , trait_plot],
                    ylim  = ylim,
                    ylab  = ylab,
                    xaxt  = 'n',
                    fsize = fsize,
                    lwd   = lwd,
                    spread.labels = spread,
                    colors = c("black" = background_col, "orange" = "orange", "blue" = "blue"),
                    main   = NULL,
                    cex.main = cex.main)
title(main = "Brownian Motion\nFast Evolutionary Rate", cex.main = cex.main)

matplot(spec_bmf, col = spec_colvec, type = "l", lty = 1, ylim = ylim_spec,
        lwd = lwd, ylab = "reflectance", xlab = "wavelength (nm)")

lines(spec_bmf[ , picks[1]], col = spec_colvec[picks[1]])
lines(spec_bmf[ , picks[2]], col = spec_colvec[picks[2]])


####################
phytools::phenogram(tree  = tree_poly,
                    x     = data_bms[ , trait_plot],
                    ylim  = ylim,
                    ylab  = ylab,
                    xaxt  = 'n',
                    fsize = fsize,
                    lwd   = lwd,
                    spread.labels = spread,
                    colors = c("black" = background_col, "orange" = "orange", "blue" = "blue"),
                    main   = NULL,
                    cex.main = cex.main)
title(main = "Brownian Motion\nSlow evolutionary Rate", cex.main = cex.main)

matplot(spec_bms, col = spec_colvec, type = "l", lty = 1, ylim = ylim_spec,
        lwd = lwd, ylab = "reflectance", xlab = "wavelength (nm)")

lines(spec_bms[ , picks[1]], col = spec_colvec[picks[1]])
lines(spec_bms[ , picks[2]], col = spec_colvec[picks[2]])

####################
phytools::phenogram(tree  = tree_poly,
                    x     = data_ouf[ , trait_plot],
                    ylim  = ylim,
                    ylab  = ylab,
                    xaxt  = 'n',
                    fsize = fsize,
                    lwd   = lwd,
                    spread.labels = spread,
                    colors = c("black" = background_col, "orange" = "orange", "blue" = "blue"),
                    main   = NULL,
                    cex.main = cex.main)
title(main = "Ornstein-Uhlenbeck\nSingle optimum", cex.main = cex.main)

matplot(spec_ouf, col = spec_colvec, type = "l", lty = 1, ylim = ylim_spec,
        lwd = lwd, ylab = "reflectance", xlab = "wavelength (nm)")
lines(spec_ouf[ , picks[1]], col = spec_colvec[picks[1]])
lines(spec_ouf[ , picks[2]], col = spec_colvec[picks[2]])


dev.off()
