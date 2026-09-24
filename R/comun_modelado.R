# Shared helpers for scripts 19 and 20: data, formulas, metrics and GBM fit.
suppressPackageStartupMessages({library(data.table); library(xgboost)})
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
N_FONDO <- 10000
ACC   <- c("dist_vias_km","dist_poblados_km","dist_rios_km")
GEO   <- c("era","amb","log_edad","edad_na","woe_uni")
TOPO  <- c("elevacion","pendiente","rugosidad")
EXPO  <- "exposicion"                          # rock exposure, script 02c
form  <- function(v) as.formula(paste("~", paste(v, collapse=" + "), "- 1"))
F_GEO <- form(GEO); F_GT0 <- form(c(GEO, TOPO)); F_GT <- form(c(GEO, TOPO, EXPO))
F_ALL <- form(c(GEO, TOPO, EXPO, ACC))

auc <- function(y,p){ r<-rank(p); n1<-sum(y==1); n0<-sum(y==0); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
tss <- function(y,p){ u <- unique(quantile(p, seq(0.5,0.999,0.001), na.rm=TRUE))
  max(sapply(u, function(t) mean(p[y==1]>=t)+mean(p[y==0]<t)-1)) }
pr_auc <- function(y,p){ o<-order(-p); y<-y[o]; tp<-cumsum(y==1); fp<-cumsum(y==0)
  sum(diff(c(0, tp/sum(y==1))) * tp/(tp+fp)) }
boyce <- function(y, p, nclas = 10) {
  br <- unique(quantile(p, seq(0,1,length.out=nclas+1), na.rm=TRUE))
  if (length(br) < 3) return(NA_real_)
  br[1] <- -Inf; br[length(br)] <- Inf
  cl <- cut(p, br, labels=FALSE); k <- length(br)-1
  pe <- sapply(1:k, function(i){ E <- sum(cl==i)/length(cl)
                                 if (E==0) NA else (sum(y[cl==i]==1)/sum(y==1))/E })
  if (sum(!is.na(pe)) < 3) return(NA_real_)
  suppressWarnings(cor(1:k, pe, method="spearman", use="complete.obs"))
}
woe_map <- function(tr){ N<-nrow(tr); Dn<-sum(tr$pres)
  w <- tr[, .(n=.N,d=sum(pres)), by=uni]; w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]
  setNames(w$Wmas,w$uni) }
apl <- function(m,u){ v<-m[u]; fifelse(is.na(v),0,as.numeric(v)) }

# Split by a logical vector (TRUE = test); WoE estimated on training only
dividir <- function(es_test){
  tr <- D[!es_test][!(pres==1 & unidad_ignea)]; te <- D[es_test]
  wm <- woe_map(tr); tr[, woe_uni := apl(wm,uni)]; te[, woe_uni := apl(wm,uni)]
  list(tr=tr, te=te, ent=rbind(tr[pres==1], tr[pres==0][sample(.N,min(.N,N_FONDO))])) }
gbm <- function(f, ent){ y <- ent$pres
  xgb.train(params=list(objective="binary:logistic",max_depth=5,eta=0.05,subsample=0.8,
                        colsample_bytree=0.8,scale_pos_weight=sum(y==0)/sum(y==1),nthread=8),
            data=xgb.DMatrix(model.matrix(f,ent),label=y), nrounds=300, verbose=0) }
pred <- function(m, f, dt) predict(m, xgb.DMatrix(model.matrix(f, dt)))
neutralizar <- function(te, tr){ x <- copy(te)
  for (v in ACC) x[[v]] <- median(tr[[v]][tr$pres==1], na.rm=TRUE); x }
# Geological baseline as a logistic GLM on a fixed design matrix, so that
# classes absent from training get a zero coefficient instead of an error
base_glm <- function(tr, te){
  Xtr <- model.matrix(~ era + amb + woe_uni - 1, tr); Xte <- model.matrix(~ era + amb + woe_uni - 1, te)
  cf <- suppressWarnings(glm.fit(Xtr, tr$pres, family=binomial()))$coefficients
  cf[is.na(cf)] <- 0; as.numeric(plogis(Xte %*% cf)) }
metricas <- function(y, p, dvias, pref) setNames(
  list(auc(y,p), tss(y,p), pr_auc(y,p), boyce(y,p), cor(p, dvias, method="spearman")),
  paste0(pref, c("_AUC","_TSS","_PR","_Boyce","_rhoVias")))
desplazamiento <- function(a, b){ ta <- a >= quantile(a,0.9); tb <- b >= quantile(b,0.9)
  100*(1 - sum(ta & tb)/sum(tb)) }
