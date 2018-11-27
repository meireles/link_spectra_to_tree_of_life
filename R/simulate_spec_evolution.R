library("phytools")
library("mvMORPH")
library("PEcAnRTM")

################################################################################
# Simulate Tree
################################################################################

# This seed outputs a reasonable tree shape
set.seed(11)

# Base tree
tree_base = phytools::pbtree(n = 50, scale = 1)

# Unset seed
set.seed(Sys.time())

# And another tree with a convergent `regime`
tree_poly = phytools::paintSubTree(tree_base, node = 90, stem = TRUE,
                                   state = "understory", anc.state = "sun")

tree_poly = phytools::paintSubTree(tree_poly, node = 61, stem = TRUE,
                                   state = "understory", anc.state = "sun")

# Pick regime colors
regime_cols = setNames( c("orange", "darkblue"),
                        c("sun", "understory") )

regime_look = setNames(rep("sun", Ntip(tree_poly)), tree_poly$tip.label)

understory = c(getDescendants(tree_poly, 90), getDescendants(tree_poly, 61))
understory = understory[understory <= Ntip(tree_poly)]

regime_look[understory] = "understory"

################################################################################
# Simulate Traits
################################################################################

prospect_mu = c("N"   = 1.6,     # structural param. N of layers
                "Cab" = 20.0,    # chlorophyll a + b in ug cm^-2
                "Car" = 12.0,    # total carotenoids in ug cm^-2
                "Cw"  = 0.02,    # equivalent water thickness (EWT) in cm
                "Cm"  = 0.01)    # LMA in g cm^-2

sigma = c("N" = 0.1, "Cab" = 1, "Car" = 0.5, "Cw"  = 0.001, "Cm"  = 0.0005)

prospect_mu_alt  = prospect_mu
prospect_mu_alt["Cab"] = 60.0

# Sigma^2
oum_sigma_sq = sigma * sigma

# Alpha
alpha       = log(2)
alpha_small = 1e-9
oum_alpha   = matrix(data = diag( c(alpha_small, alpha, alpha_small, alpha_small, alpha_small)),
                     nrow = 5, ncol = 5,
                     dimnames = list(names(prospect_mu),
                                     names(prospect_mu)) )
# Theta
oum_theta = rbind("sun"        = prospect_mu,
                  "understory" = prospect_mu_alt)

example_params_oum = list(ntraits      = 5,
                          sigma        = diag(oum_sigma_sq),
                          alpha        = oum_alpha,
                          theta        = oum_theta,
                          names_traits = names(prospect_mu),
                          vcv          = "fixedRoot")

example_data_oum     = mvMORPH::mvSIM(tree  = tree_poly,
                                      param = example_params_oum,
                                      model = "OUM",
                                      nsim  = 1)

################################################################################
# Simulate Spectra
################################################################################

spec = apply(example_data_oum, 1, function(x){
    PEcAnRTM::prospect(x, version = "5")[ , 1]
})
spec_anc = PEcAnRTM::prospect(prospect_mu, version = "5")[ , 1]

rownames(spec)     = 400:2500
rownames(spec_anc) = 400:2500

spec     = spec[rownames(spec) %in% seq(400, 2500, 10), ]
spec_anc = spec_anc[rownames(spec_anc) %in% seq(400, 2500, 10), ]


################################################################################
# Estimate model prob
################################################################################

# Looping like this is kinda stupid.
# multivar takes too long. switch packages...

params = list("optimization" = "subplex")
models_aicw = sapply(seq(nrow(spec)), function(x){

    bm1 = mvMORPH::mvBM(tree_poly, data = t(spec)[ , x], model = "BM1", param = params)
    bmm = mvMORPH::mvBM(tree_poly, data = t(spec)[ , x], model = "BMM", param = params)
    ou1 = mvMORPH::mvOU(tree_poly, data = t(spec)[ , x], model = "OU1", param = params)
    oum = mvMORPH::mvOU(tree_poly, data = t(spec)[ , x], model = "OUM", param = params)

    aicw(list(bm1, bmm, ou1, oum))$aicweights
})


################################################################################
# Plot
################################################################################

png("figures/ou_simulation_and_inference.png", width = 6, height = 6, units = "in", res = 800)

mar = c(4, 4, 3, 2)
par(mfrow = c(2, 2), mar = mar)

########################################
# Tree
########################################
phytools::plotSimmap(tree_poly, colors = regime_cols, lwd = 1.5, fsize = 1e-9,
                     mar = mar)

legend("bottomleft", names(regime_cols), fill = regime_cols, bty = "n")

########################################
# Traitgram
########################################

phenogram(tree = tree_poly, x = example_data_oum[ , "Cab"],fsize = 1e-9, spread.labels = F,
          colors = regime_cols,lwd = 0.8, ylab = "Chlorophyll content", xlab = "time")

########################################
# Spectra
########################################

## here we can choose to plot the spectra themselves or their difference from the
## ancestral state

## Spectra
#w_range = 1:81
# matplot(y = spec[ w_range, ], x = rownames(spec)[w_range],type = "l", col = regime_cols[regime_look],
#         ylab = "reflectance", xlab = "wavelengths", lty = 1, lwd = 0.5)

## Difference
spec_diff = spec - as.vector(spec_anc)
matplot(y = spec_diff, x = rownames(spec_diff),type = "l", col = regime_cols[regime_look],
        ylab = "Diff. from ancestral spectrum", xlab = "wavelengths", lty = 1, lwd = 0.5)
abline(h = 0, lty = 2)

########################################
## Model weights
########################################

model_names = setNames(c("purple", "green", "orange", "red"),
                       c("BM1", "BMM", "OU1", "OUM"))

plot(y = models_aicw[4 , ], x = rownames(spec), type = "l", lty = 1, lwd = 1.5,
        col = "black", ylab = "AIC weight OUM model", xlab = "wavelengths")

dev.off()
