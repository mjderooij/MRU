setwd("~/surfdrive/LogitMDA/mru/Analysis")
library(lmap)
library(nnet)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(GGally, quietly = TRUE, warn.conflicts = FALSE)
library(mvtnorm)
library(ggpubr)
library(metR)
#library(tidyverse)
library(RColorBrewer)
library(jtools)

# DATA
load("~/surfdrive/LogitMDA/dpes02.Rdata")
mydat3= mydat2[, 14:19]
mydat3 = mydat3[complete.cases(mydat3), ]
X = mydat3[ , 1:5]; y = mydat3[, 6]
G = class.ind(y)
G = G[, c(1,2,3,4,5,7,9,10)]
X = X[which(rowSums(G) != 0), ]
X = as.matrix(X)
#X = X - outer(rep(1,nrow(X)), c(4,4,4,4,6))
X = scale(X)
G = G[which(rowSums(G) != 0), ]
colnames(G) = c("PvdA", "CDA", "VVD", "D66", "GL" , "CU", "LPF", "SP")
y = apply(G, 1, which.max)
y = factor(y, 
           levels = c("1","2" ,"3", "4", "5", "6", "7", "8"), 
           labels = c("PvdA", "CDA", "VVD", "D66", "GL" , "CU", "LPF", "SP"))

xnames = colnames(X)
ynames = colnames(G)

# ANALYSES
set.seed(1234)
R = 200
M1results = vector(mode = "list", length = R)
M2results = vector(mode = "list", length = R)
M3results = vector(mode = "list", length = R)
M4results = vector(mode = "list", length = R)
M5results = vector(mode = "list", length = R)

M1results[[1]] = mru(y = y, X = X, S = 1, start = "da")
M2results[[1]] = mru(y = y, X = X, S = 2, start = "da")
M3results[[1]] = mru(y = y, X = X, S = 3, start = "da")
M4results[[1]] = mru(y = y, X = X, S = 4, start = "da")
M5results[[1]] = mru(y = y, X = X, S = 5, start = "da")

for(r in 2:R){
  M1results[[r]] = mru(y = y, X = X, S = 1, start = "random")
  M2results[[r]] = mru(y = y, X = X, S = 2, start = "random")
  M3results[[r]] = mru(y = y, X = X, S = 3, start = "random")
  M4results[[r]] = mru(y = y, X = X, S = 4, start = "random")
  M5results[[r]] = mru(y = y, X = X, S = 5, start = "random")
}
#save(M1results, M2results, M3results, M4results, M5results, file = "dpes02results.Rdata")

library(ggplot2)
dev1 = dev2 = dev3 = dev4 = dev5 = rep(NA, 100)
for(r in 1:R){
  dev1[r] = M1results[[r]]$deviance
  dev2[r] = M2results[[r]]$deviance
  dev3[r] = M3results[[r]]$deviance
  dev4[r] = M4results[[r]]$deviance
  dev5[r] = M5results[[r]]$deviance
}
df = data.frame(repl = rep(1:R, 5), dimensionality = as.factor(rep(1:5, each = R)), deviance = c(dev1, dev2, dev3, dev4, dev5))
p = ggplot(df, aes(dimensionality, deviance)) + geom_boxplot() + theme_apa()
p
#ggsave("~/surfdrive/multldm/mrmdu/paper/figures/dpes02results.pdf", plot = p)

# DIMENSIONALITY SELECTION
out.dpes1 = M1results[[which.min(dev1)]]
out.dpes2 = M2results[[which.min(dev2)]]
out.dpes3 = M3results[[which.min(dev3)]]
out.dpes4 = M4results[[which.min(dev4)]]
out.dpes5 = M5results[[which.min(dev5)]]

fit = matrix(NA, 5, 4)
fit[ , 1] = c(1,2,3, 4, 5)
fit[1, 2] = out.dpes1$deviance
fit[2, 2] = out.dpes2$deviance
fit[3, 2] = out.dpes3$deviance
fit[4, 2] = out.dpes4$deviance
fit[5, 2] = out.dpes5$deviance
fit[1, 3] = ncol(X) * 1 + ncol(G) * 1 
fit[2, 3] = ncol(X) * 2 + ncol(G) * 2 - 2*1/2 
fit[3, 3] = ncol(X) * 3 + ncol(G) * 3 - 3*2/2
fit[4, 3] = ncol(X) * 4 + ncol(G) * 4 - 4*3/2
fit[5, 3] = ncol(X) * 5 + ncol(G) * 5 - 5*4/2
fit[, 4] = fit[, 2] + 2* fit[, 3]
colnames(fit) = c("Dimensionality", "Deviance", "#params", "AIC")
fit

# VARIABLE SELECTION
BB = out.dpes2$B
VV = out.dpes2$V

out.dpes2.1 = fastmru(G = G, X = X[, -1], B = BB[-1, ], V = VV, DCRIT = 1e-05)
out.dpes2.2 = fastmru(G = G, X = X[, -2], B = BB[-2, ], V = VV, DCRIT = 1e-05)
out.dpes2.3 = fastmru(G = G, X = X[, -3], B = BB[-3, ], V = VV, DCRIT = 1e-05)
out.dpes2.4 = fastmru(G = G, X = X[, -4], B = BB[-4, ], V = VV, DCRIT = 1e-05)
out.dpes2.5 = fastmru(G = G, X = X[, -5], B = BB[-5, ], V = VV, DCRIT = 1e-05)

fit2 = matrix(NA, 6, 6)
fit2[ , 1] = c(0, 1, 2, 3, 4, 5)
fit2[1, 2] = out.dpes2$deviance
fit2[2, 2] = out.dpes2.1$deviance
fit2[3, 2] = out.dpes2.2$deviance
fit2[4, 2] = out.dpes2.3$deviance
fit2[5, 2] = out.dpes2.4$deviance
fit2[6, 2] = out.dpes2.5$deviance
fit2[1 , 3] = (ncol(X)) * 2 + ncol(G) * 2 - 2*1/2 
fit2[2:6 , 3] = (ncol(X) -1) * 2 + ncol(G) * 2 - 2*1/2 
fit2[ , 4] = fit2[, 2] + 2* fit2[, 3]
fit2[, 5] = fit2[, 2] - fit2[1, 2]
fit2[, 6] = round(1 - pchisq(fit2[, 5], df = 2), digits = 3)
colnames(fit2) = c("Left out X", "Deviance", "#params", "AIC", "LRT", "p")
rownames(fit2) = c("-", colnames(X)) 
fit2

save(X, y, G, out.dpes2, file = "dpes_finresult.Rdata")

# VISUALIZATION
p1 = plot(out.dpes2, ocol = y)
ggsave("~/surfdrive/LogitMDA/mru/Paper/figures/dpes.pdf", plot = p1, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

# CLASSIFICATION PERFORMANCE
U = out.dpes2$U; V = out.dpes2$V
D = outer(diag(U %*% t(U)), rep(1, 8)) + outer(rep(1,275), diag(V %*% t(V))) - 2* U %*% t(V)
D = sqrt(D)
PR = exp(-D)/rowSums(exp(-D))
yhat = apply(PR, 1, which.max)
tab1 = table(y, yhat)
tab1

100*sum(diag(tab1))/sum(tab1)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# squared distance model 
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
ipc.deviance <- function(pars, G, X){
  # Deviance of IPC model for LIVER data
  # pars: initial values of parameters
  # G: indicator matrix representing the response classes(n  by J)
  # X: predictor variables (n by p + 1); should include a vector of ones
  # CALL: stats <- optim(pars, ipc.deviance2, NULL, Y, X, Z, method ="BFGS", 
  #                      control = list(trace = 2, maxit = 100), hessian = F)
  # Change to hessian = T if se's need to be computed.
  # \copyright M. de Rooij, 23-5-2013                 \
  # ------------------------------------------------------------------
  # extract matrix B from the pars
  I = nrow(G)
  C = ncol(G)
  p = ncol(X)
  B = matrix(pars[1:(2*p)], p, 2)
  # create Z - matrix with coordinates of class points
  # TRANSLATION, SCALING, AND ROTATION
  V = matrix(0, C, 2)
  #V[2,1] = V[3,2] = 1
  V[2:C,1] = pars[(2*p+1):(2*p+7)]
  V[3:C,2] = pars[(2*p+8):(2*p+13)]
  
  # Make the linear predictors and compute the squred distances
  N = X %*% B
  
  v2 = rowSums(V ^ 2)
  D =  - 2 * N %*% t(V) + rep(1, I) %o% v2  
  
  # Compute probabilities and deviance
  P = exp(-D)
  sp = rowSums(P)
  P = (1 / sp) * P
  dev = -2 * sum(G * log(P))
}

XX = cbind(1,X)
# XX = cbind(1,scale(X))

pars0 = c(0, out.dpes2$B[,1], 0, out.dpes2$B[,2], out.dpes2$V[2:8, 1], out.dpes2$V[3:8, 2])
out = optim(pars0,ipc.deviance, G = G, X = XX, method = "BFGS", control = list(trace = 2, maxit = 1000), hessian = F)
dev2.2 = rep(NA, 100); dev2.2[1] = out$value
set.seed(1234)
for(r in 2:100){
  pars0 = rnorm(length(pars0))
  out2 = optim(pars0,ipc.deviance, G = G, X = XX, method = "BFGS", control = list(trace = 0, maxit = 1000), hessian = F)
  dev2.2[r] = out2$value
  if(out2$value < out$value){out = out2}
}

# AIC
out$value + 2* length(out$par)

P = ncol(XX); C = ncol(G)
B = matrix(out$par[1:(2*P)], P, 2)
U = X %*% B
V = matrix(0, C, 2)
V[2:C,1] = out$par[(2*P+1):(2*P+7)]
V[3:C,2] = out$par[(2*P+8):(2*P+13)]

# percentage correct
D2 = outer(diag(U %*% t(U)), rep(1, 8)) + outer(rep(1,275), diag(V %*% t(V))) - 2* U %*% t(V)
PR2 = exp(-D2)/rowSums(exp(-D2))
yhat2 = apply(PR2, 1, which.max)
tab2 = table(y, yhat2)
100*sum(diag(tab2))/sum(tab2)

# plot
VV = V - outer(rep(1,8), B[1, ])
# rotation to principal axes
UU = XX[, -1] %*% B[-1, ]
R = eigen(t(UU)%*%UU)$vectors
B = B[-1, ] %*% R
UU = UU %*% R
VV = VV %*% R

NN = as.data.frame(UU)
colnames(NN) = c("Dim1", "Dim2")
VV = as.data.frame(VV)  
colnames(VV) = c("Dim1", "Dim2")
P = nrow(B)
Xo = XX[, -1]
# for solid line
MCx1 <- data.frame(labs=character(),
                   varx = integer(),
                   Dim1 = double(),
                   Dim2 = double(), stringsAsFactors=FALSE)
# for markers
MCx2 <- data.frame(labs=character(),
                   varx = integer(),
                   Dim1 = double(),
                   Dim2 = double(), stringsAsFactors=FALSE)

ll = 0
lll = 0
for(pp in 1:P){
  b = matrix(B[pp , ], 2, 1)
  # solid line
  minx = min(Xo[, pp])
  maxx = max(Xo[, pp])
  m.x1 = c(minx,maxx)
  markers1 = matrix(m.x1, 2, 1) 
  markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
  MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), pp)
  MCx1[(ll + 1): (ll + 2), 2] = pp
  MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1 
  ll = ll + 2
  # markers
  m.x2 = pretty(Xo[, pp]) 
  m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
  l.m = length(m.x2)
  markers2 = matrix(m.x2, l.m, 1) 
  markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
  MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2) 
  MCx2[(lll + 1): (lll + l.m), 2] = pp
  MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2 
  lll = lll + l.m
} # loop p

p2 = ggplot() +
  geom_point(data = VV, aes(x = Dim1, y = Dim2), colour = "darkgreen", size = 3) +
  geom_point(data = NN, aes(x = Dim1, y = Dim2), colour = "grey", size = 1, show.legend = FALSE) + 
  xlab("Dimension 1") + 
  ylab("Dimension 2") 

p2 = p2 + geom_text_repel(data = VV, aes(x = Dim1, y = Dim2), label = colnames(G))

xcol = "lightskyblue"
p2 = p2 + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = xcol, linetype = 3) +
  geom_line(data = MCx1, aes(x = Dim1, y = Dim2, group = varx), col = xcol, linewidth = 1) + 
  geom_point(data = MCx2, aes(x = Dim1, y = Dim2), col = xcol) + 
  geom_text(data = MCx2, aes(x = Dim1, y = Dim2, label = labs), nudge_y = -0.08, size = 1.5)

a = (max(abs(c(ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range, ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)))) + 0.1

idx1 = apply(abs(B), 1, which.max)
t = s = rep(NA,(P))
for(pp in 1:(P)){
  t[(pp)] = (a *1.1)/(abs(B[pp,idx1[(pp)]])) * B[pp,-idx1[(pp)]]
  s[(pp)] = sign(B[pp,idx1[(pp)]])
}
CC = cbind(idx1, t, s)
bottom = which(CC[, "idx1"] == 2 & CC[, "s"] == -1)
top =  which(CC[, "idx1"] == 2 & CC[, "s"] == 1)
right = which(CC[, "idx1"] == 1 & CC[, "s"] == 1)
left = which(CC[, "idx1"] == 1 & CC[, "s"] == -1)

p2 = p2 + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = xnames[bottom], 
                             sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = xnames[top]))
p2 = p2 + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = xnames[left], 
                             sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = xnames[right]))  
p2 = p2 + coord_fixed()
p2 = p2 + theme_bw() 
p2 = p2 + theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) 


p2
ggsave("~/surfdrive/logitMDA/mru/paper/figures/dpesplot2.pdf", plot = p2, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)


# multinomial logistic regression

df = as.data.frame(cbind(y, X))
mlr.out = multinom(y ~ E + ID + AS + C + LR, data = df, trace = FALSE)
