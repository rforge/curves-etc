library("plugdensity")

data(geyser)
str(pd.geys <- plugin.density(geyser$waiting))
