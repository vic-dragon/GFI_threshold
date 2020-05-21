setwd("/home/syhwang/Simulation7")
library(segmented)
library(GA)
source("GenerateSamples.R")
source("ModelSelection.R")
source("Densities.R")
source("Fiducial.R")
source("pvaluemodel.R")
source("Simulation.R")

n=20
m=50
sigma=0.1

store.f1<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=1,sigma=sigma)
  store.f1[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f1_finish')

store.f2<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=2,sigma=sigma)
  store.f2[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f2_finish')

store.f3<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=3,sigma=sigma)
  store.f3[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f3_finish')

img_name=paste0("n",n,"m",m,"s0",sigma*10,'.Rdata')
save.image(file=img_name)

sigma=0.3

store.f1<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=1,sigma=sigma)
  store.f1[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f1_finish')

store.f2<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=2,sigma=sigma)
  store.f2[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f2_finish')

store.f3<-NULL
for(i in 1:200){
  dataset = gen.sample(n=n,m=m,fnct=3,sigma=sigma)
  store.f3[[i]] = Fid.spline(dataset, old.threshold = 0.5, sd.threshold = 0.1,burn.in=1000, s.size=2000)
}
print('f3_finish')


img_name=paste0("n",n,"m",m,"s0",sigma*10,'.Rdata')
save.image(file=img_name)
