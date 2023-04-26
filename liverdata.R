library(nnet)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(GGally, quietly = TRUE, warn.conflicts = FALSE)
library(mvtnorm)
library(lmap)
library(ggpubr)
library(metR)
library(tidyverse)
library(RColorBrewer)
library(jtools)

load("~/surfdrive/MultinomialBook/Datasets/liverdata/liverdata.RData")
colnames(liverdat) = c("Disease", "AS", "AL", "GD", "X?")
X = liverdat[ ,-c(1,5)]
xnames= colnames(X)
X.df = data.frame(X)
y = liverdat[ , 1]
y = factor(y, levels = c("1","2" ,"3", "4"), labels = c("AVH","PCH","ACH","PNC"))
G = class.ind(y)
# colnames(G) = c("AVH","PCH","ACH","PNC")
# G = G[, c(4,1,2,3)]
XX = log(X)
X = scale(XX, center = TRUE, scale = FALSE)

# analysis with distacne model
set.seed(1234)

M1results = vector(mode = "list", length = 100)
M2results = vector(mode = "list", length = 100)
M3results = vector(mode = "list", length = 100)

M1results[[1]] = mru(y = y, X = X, S = 1, start = "da")
M2results[[1]] = mru(y = y, X = X, S = 2, start = "da")
M3results[[1]] = mru(y = y, X = X, S = 3, start = "da")


for(r in 2:100){
  M1results[[r]] = mru(y = y, X = X, S = 1, start = "random")
  M2results[[r]] = mru(y = y, X = X, S = 2, start = "random")
  M3results[[r]] = mru(y = y, X = X, S = 3, start = "random")
}

dev1 = dev2 = dev3 = rep(NA, 100)
for(r in 1:100){
  dev1[r] = M1results[[r]]$deviance
  dev2[r] = M2results[[r]]$deviance
  dev3[r] = M3results[[r]]$deviance
}

df = data.frame(repl = rep(1:100, 3), dimensionality = as.factor(rep(1:3, each = 100)), deviance = c(dev1, dev2, dev3))
ggplot(df, aes(dimensionality, deviance)) + geom_boxplot() + theme_bw()

mru.out1 = M1results[[which.min(dev1)]]
mru.out2 = M2results[[which.min(dev2)]]
mru.out3 = M3results[[which.min(dev3)]]

source("~/surfdrive/LogitMDA/lmap-package/new/R/plot.mru.R")
plt.mru = plot(mru.out2, ocol = NA)

# color objects by residuals
source("~/surfdrive/LogitMDA/mru/Analysis/fitted.mru.R")
Ghat1 = fitted.mru(mru.out2)$Ghat
resid.mru = 2* rowSums(mru.out2$G * log(1/Ghat1))
df.mru = cbind.data.frame(mru.out2$U, resid.mru); colnames(df.mru) = c("u1","u2","residual")
plt.mru = plt.mru + geom_point(data = df.mru, aes(x = u1, y = u2, colour = residual), show.legend = FALSE) + 
  scale_color_continuous(guide = "legend",limits = c(0,9), low="grey", high="red") 
plt.mru
ggsave("~/surfdrive/LogitMDA/mru/Paper/figures/liver1.pdf", plot = plt.mru, width = 8.3) #, width = 8.3, height = 8.3, units = "in", limitsize = FALSE)

# color by output class
df.mru = cbind.data.frame(mru.out2$U, mru.out2$y); colnames(df.mru) = c("u1","u2","class")
plt.mru = plt.mru + geom_point(data = df.mru, aes(x = u1, y = u2, colour = class), show.legend = FALSE) + 
  scale_color_manual(values = c("AVH" = "green", "PCH" ="orange", "ACH"="steelblue", "PNC" = "red")) 
plt.mru



#+ 
  # scale_color_continuous(type = "viridis") #+ 
  #scale_size_continuous(limits = c(0,9)) #+ labs(title = "distance model")

# analysis with squared distance model 
source("~/surfdrive/LogitMDA/mru/Analysis/ipc4.R")

ipc.out = ipc4(G, X)
dev2.2 = rep(NA, 100); dev2.2[1] = ipc.out$deviance
set.seed(1234)
for(r in 2:100){ # niet nodig - elke keer zelfde oplossing
  out2 = ipc4(G, X, random = TRUE)
  dev2.2[r] = out2$deviance
  if(out2$deviance < ipc.out$deviance){ipc.out = out2}
}

R = eigen(t(ipc.out$U)%*%ipc.out$U)$vectors
ipc.out$U = ipc.out$U %*% R
ipc.out$B = ipc.out$B %*% R
ipc.out$V = ipc.out$V %*% R

plt.ipc = plot(ipc.out, ocol = NA)

# color by residuals
Ghat2 = fitted.mru(ipc.out, squared = TRUE)$Ghat
resid.ipc = 2* rowSums(ipc.out$G * log(1/Ghat2))
df.ipc = cbind.data.frame(ipc.out$U, resid.ipc); colnames(df.ipc) = c("u1","u2","residual")
plt.ipc = plt.ipc + geom_point(data = df.ipc, aes(x = u1, y = u2, colour = residual), show.legend = FALSE) + 
  scale_color_continuous(guide = "legend",limits = c(0,9), low="grey", high="red") #+ 
  # scale_size_continuous(limits = c(0,9)) #+ labs(title = "squared distance model")
plt.ipc
# ggsave("~/surfdrive/LogitMDA/mru/Paper/figures/liver2.pdf", plot = plt.ipc, width = 8.3) #, width = 8.3, height = 8.3, units = "in", limitsize = FALSE)

# color by output class
df.ipc = cbind.data.frame(ipc.out$U, mru.out2$y); colnames(df.ipc) = c("u1","u2","class")
plt.ipc = plt.ipc + geom_point(data = df.ipc, aes(x = u1, y = u2, colour = class), show.legend = FALSE) + 
  scale_color_manual(values = c("AVH" = "green", "PCH" ="orange", "ACH"="steelblue", "PNC" = "red")) 
plt.ipc = plt.ipc + geom_point(aes(x = 4.295, y = -3.097), colour = "red", size = 4, shape = "diamond")

rc = -1 / ((ipc.out$V[4, 2] - ipc.out$V[3, 2]) / (ipc.out$V[4, 1] - ipc.out$V[3, 1]))
xy = ipc.out$U[180, ]
plt.ipc = plt.ipc + geom_abline(intercept = (xy[2] - rc * xy[1]), slope = rc, col = "grey", linetype = 3)
plt.ipc

plt.ipc = plt.ipc + geom_point(aes(x = 6.53989143, y = 0.438196894), colour = "red", size = 4, shape = "diamond")

plt2 = ggpubr::ggarrange(plt.mru, plt.ipc, nrow = 2)
ggsave("~/surfdrive/LogitMDA/mru/Paper/figures/liverclass.pdf", plot = plt2, width = 8.3) #, width = 8.3, height = 8.3, units = "in", limitsize = FALSE)

# P = ncol(X); C = ncol(G)
# B = matrix(out$par[1:(2*P)], P, 2)
# U = X %*% B
# V = matrix(0, C, 2)
# V[2,1] = V[3,2] = 1
# V[3:C,1] = out$par[(2*P+1):(2*P+2)]
# V[4:C,2] = out$par[(2*P+3)]
# 
# # plot
# VV = V - outer(rep(1,4), B[1, ])
# UU = X[, -1] %*% B[-1, ]
# R = eigen(t(UU)%*%UU)$vectors
# B = B[-1, ] %*% R
# UU = UU %*% R
# VV = VV %*% R
# 
# NN = as.data.frame(UU)
# colnames(NN) = c("Dim1", "Dim2")
# 
# VV = as.data.frame(V)  
# colnames(VV) = c("Dim1", "Dim2")
# 
# P = nrow(B)
# Xo = X[ , -1]
# 
# # for solid line
# MCx1 <- data.frame(labs=character(),
#                    varx = integer(),
#                    Dim1 = double(),
#                    Dim2 = double(), stringsAsFactors=FALSE)
# # for markers
# MCx2 <- data.frame(labs=character(),
#                    varx = integer(),
#                    Dim1 = double(),
#                    Dim2 = double(), stringsAsFactors=FALSE)
# 
# ll = 0
# lll = 0
# for(pp in 1:P){
#   b = matrix(B[pp , ], 2, 1)
#   # solid line
#   minx = min(Xo[, pp])
#   maxx = max(Xo[, pp])
#   m.x1 = c(minx,maxx)
#   markers1 = matrix(m.x1, 2, 1) 
#   markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
#   MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), pp)
#   MCx1[(ll + 1): (ll + 2), 2] = pp
#   MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1 
#   ll = ll + 2
#   # markers
#   m.x2 = pretty(Xo[, pp]) 
#   m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
#   l.m = length(m.x2)
#   markers2 = matrix(m.x2, l.m, 1) 
#   markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
#   MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2) 
#   MCx2[(lll + 1): (lll + l.m), 2] = pp
#   MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2 
#   lll = lll + l.m
# } # loop p
# 
# p2 = ggplot() +
#   geom_point(data = VV, aes(x = Dim1, y = Dim2), colour = "darkgreen", size = 3) +
#   geom_point(data = NN, aes(x = Dim1, y = Dim2, color = y), size = 1, show.legend = FALSE) + 
#   xlab("Dimension 1") + 
#   ylab("Dimension 2") 
# 
# p2 = p2 + geom_text(data = VV, aes(x = Dim1, y = Dim2), 
#                     label = colnames(G),
#                     vjust = 0, nudge_y = -0.5)
# 
# 
# xcol = "lightskyblue"
# p2 = p2 + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = xcol, linetype = 3) +
#     geom_line(data = MCx1, aes(x = Dim1, y = Dim2, group = varx), col = xcol, linewidth = 1) + 
#     geom_point(data = MCx2, aes(x = Dim1, y = Dim2), col = xcol) + 
#     geom_text(data = MCx2, aes(x = Dim1, y = Dim2, label = labs), nudge_y = -0.08, size = 1.5)
#   
# a = (max(abs(c(ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range, ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))) + 0.1
#   
# idx1 = apply(abs(B), 1, which.max)
# t = s = rep(NA,(P))
# for(pp in 1:(P)){
#   t[(pp)] = (a * 1.1)/(abs(B[pp,idx1[(pp)]])) * B[pp,-idx1[(pp)]]
#   s[(pp)] = sign(B[pp,idx1[(pp)]])
# }
# CC = cbind(idx1, t, s)
# bottom = which(CC[, "idx1"] == 2 & CC[, "s"] == -1)
# top =  which(CC[, "idx1"] == 2 & CC[, "s"] == 1)
# right = which(CC[, "idx1"] == 1 & CC[, "s"] == 1)
# left = which(CC[, "idx1"] == 1 & CC[, "s"] == -1)
# 
# 
# p2 = p2 + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = xnames[bottom], 
#                              sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = xnames[top]))
# p2 = p2 + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = xnames[left], 
#                              sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = xnames[right]))  
# p2 = p2 + coord_fixed()
# 
# p2 = p2 + theme_bw()  
# p2 = p2 + theme(axis.line = element_line(colour = "black"),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.border = element_blank(),
#                 panel.background = element_blank()) 
# 
# p2
