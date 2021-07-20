X = expand.grid(a = seq(1e-1,5000,length.out = 20), b = seq(1e-1,5000,length.out = 20))
head(X)
library(plotly)
dat = data.frame(x = X[,1], y = X[,2], z = apply(X, 1, function(x) beta(x[1],x[2])))

fig = plot_ly(dat, x = ~x, y = ~y, z = ~z, 
              type = 'scatter3d', mode = 'markers', 
              marker = list(size = 2.))
fig
orca(fig, "figs/R_beta_issue.png", scale = 2)
