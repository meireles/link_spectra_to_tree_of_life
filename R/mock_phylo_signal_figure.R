library("phytools")
library("geiger")

set.seed(8765)
phy      = ape::rcoal(9)
phy      = ape::ladderize(phy, right = FALSE)
ord      = phy$tip.label[phy$edge[ phy$edge[ , 2] <= Ntip(phy), 2]]

l_trans  = geiger::rescale(phy, "lambda")
phy_l_70 = l_trans(0.7)
phy_l_30 = l_trans(0.3)
phy_l_00 = l_trans(0.0)
phylist  = list(phy = phy, phy_l_70 = phy_l_70,
                phy_l_30 = phy_l_30, phy_l_00 = phy_l_00)

set.seed(7)
traits   = lapply(phylist, phytools::fastBM, a = 2, sig2 = 0.1, nsim = 1)
traits   = lapply(traits, `[`, ord)

pl = lapply(traits, phytools::phylosig,
            tree = phylist$phy, method = "lambda", test = FALSE)

bk = lapply(traits, phytools::phylosig,
            tree = phylist$phy, method = "K", test = FALSE)

########################################
# lambda
########################################

png("figures/lambda_trees.png", width = 7.5, height = 3.5,
    units = "in", res = 300)

par(mfrow = c(1, 4), family = "sans")

plot(phylist$phy, type = "phylogram", show.tip.label = F,
     direction = "upwards", main = "Lambda = 1")
plot(phylist$phy_l_70, type = "phylogram", show.tip.label = F,
     direction = "upwards", main = "Lambda = 0.7")
plot(phylist$phy_l_30, type = "phylogram", show.tip.label = F,
     direction = "upwards", main = "Lambda = 0.3")
plot(phylist$phy_l_00, type = "phylogram", show.tip.label = F,
     direction = "upwards", main = "Lambda = 0")

dev.off()

########################################
# K
########################################

y     = traits[c("phy", "phy_l_00")]
ylim  = range(unlist(y))
mar1  = c(1, 3, 3, 3)
mar2  = c(1, 3, 0, 3)

png("figures/phylosig_k.png", width = 7.5, height = 3.5,
    units = "in", res = 300)

par(mfcol = c(2, 2), mar = mar1, oma = c(1, 1, 1, 1), family = "sans")

plot(x = seq(Ntip(phy)),
     y = y$phy,
     axes = FALSE, ann = FALSE, pch = 16, cex = 2, ylim = ylim, xpd = TRUE)
title(paste("K =", round(bk$phy, digits = 2)), line = 2)
axis(2, line = 1, cex.axis = 0.8)
plotTree(phylist$phy, direction = "upwards", fsize = 1e-16, mar = mar2)
mtext("a", side = 1, adj = 0, outer = T, font = 2, line = -0.5)

par(mar = mar1)
plot(x = seq(Ntip(phy)),
     y = y$phy_l_00,
     axes = FALSE, ann = FALSE, pch = 16, cex = 2, , ylim = ylim, xpd = TRUE)
title(paste("K =", round(bk$phy_l_00, digits = 2)), line = 2)
axis(2, line = 1, cex.axis = 0.8)
plotTree(phylist$phy, direction = "upwards", fsize = 1e-16 ,mar = mar2)
mtext("b", side = 1, adj = 0.5, outer = T, font = 2, line = -0.5)

dev.off()
