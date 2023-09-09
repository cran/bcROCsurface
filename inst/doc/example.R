## ----check pandoc, eval=TRUE, echo=FALSE--------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
      !rmarkdown::pandoc_available("1.14")) {
    warning(call. = FALSE, "These vignettes assume rmarkdown and pandoc version 1.14.  These were not found. Older versions will not work.")
    knitr::knit_exit()
  }

## ---- eval=FALSE, install_cran------------------------------------------------
#  install.packages("bcROCsurface")

## ---- eval=TRUE, loadlib------------------------------------------------------
library(bcROCsurface)

## ---- eval=FALSE, loaddat-----------------------------------------------------
#  data(EOC)

## ---- eval=FALSE, loaddat2----------------------------------------------------
#  head(EOC)

## ----preDiseaseFULL, eval=TRUE------------------------------------------------
dise_full <- pre_data(EOC$D.full, EOC$CA125)

## ---- eval=FALSE--------------------------------------------------------------
#  head(dise_full$dise)

## ---- eval=FALSE--------------------------------------------------------------
#  head(dise_full$dise_vec)

## ----setup, echo=FALSE--------------------------------------------------------
library(knitr)
library(rgl)
knit_hooks$set(webgl = hook_webgl)

## ----ROCsFULL, webGL = TRUE---------------------------------------------------
dise_vec_full <- dise_full$dise_vec
rocs(method = "full", diag_test =  EOC$CA125, dise_vec = dise_vec_full,
     ncp = 30, ellipsoid = TRUE, cpst = c(-0.56, 2.31))

## ----vusfull, eval = FALSE----------------------------------------------------
#  vus_mar("full", diag_test = EOC$CA125, dise_vec = dise_vec_full, ci = TRUE)

## ----disease, eval = TRUE-----------------------------------------------------
dise_na <- pre_data(EOC$D, EOC$CA125)
dise_vec_na <- dise_na$dise_vec
dise_fact_na <- dise_na$dise
rho_out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC, 
                      test = TRUE)

## ----ROCsFI, webgl = TRUE-----------------------------------------------------
rocs(method = "fi", diag_test =  EOC$CA125, dise_vec = dise_vec_na, 
     veri_stat = EOC$V, rho_est = rho_out, ncp = 30, ellipsoid = TRUE, 
     cpst = c(-0.56, 2.31))

## ----vusfi,  eval=FALSE-------------------------------------------------------
#  vus_mar(method = "fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#          veri_stat = EOC$V, rho_est = rho_out, ci = TRUE)

## ----ROCsMSI, webgl = TRUE----------------------------------------------------
rocs(method = "msi", diag_test = EOC$CA125, dise_vec = dise_vec_na, 
     veri_stat = EOC$V, rho_est = rho_out, ncp = 30,
     ellipsoid = TRUE, cpst = c(-0.56, 2.31))

## ----vusmsi,  eval=FALSE------------------------------------------------------
#  vus_mar(method = "msi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#          veri_stat = EOC$V, rho_est = rho_out, ci = TRUE)

## ----verification, eval = TRUE------------------------------------------------
pi_out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, model = "logit", 
                test = TRUE, trace = TRUE)

## ----ROCsIPW, webgl = TRUE----------------------------------------------------
rocs(method = "ipw", diag_test = EOC$CA125, dise_vec = dise_vec_na, 
     veri_stat = EOC$V, pi_est = pi_out, ncp = 30,
     ellipsoid = TRUE, cpst = c(-0.56, 2.31))

## ----vusipw, eval=FALSE-------------------------------------------------------
#  vus_mar(method = "ipw", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#          veri_stat = EOC$V, pi_est = pi_out, ci = TRUE)

## ----ROCsSPE, webgl = TRUE----------------------------------------------------
rocs(method = "spe", diag_test = EOC$CA125, dise_vec = dise_vec_na, 
     veri_stat = EOC$V, rho_est = rho_out,
     pi_est = pi_out, ncp = 30, ellipsoid = TRUE, cpst = c(-0.56, 2.31))

## ----vusspe, eval=FALSE-------------------------------------------------------
#  vus_mar(method = "spe", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#          veri_stat = EOC$V, rho_est = rho_out,
#          pi_est = pi_out, ci = TRUE)

## ----1nn, eval = TRUE---------------------------------------------------------
x_mat <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
rho_1nn <- rho_knn(x_mat, dise_vec = dise_vec_na, veri_stat = EOC$V, k = 1, 
                   type = "mahala")

## ----ROCs1NN, webgl = TRUE----------------------------------------------------
rocs("knn", diag_test = EOC$CA125, dise_vec_na, veri_stat = EOC$V, 
     rho_est = rho_1nn, ncp = 30, ellipsoid = TRUE,
     cpst = c(-0.56, 2.31))

## ----vus1nn, eval = FALSE-----------------------------------------------------
#  vus_mar(method = "knn", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#          veri_stat = EOC$V, rho_est = rho_1nn, ci = TRUE,
#          parallel = TRUE)

## ----cvKnn, eval = TRUE-------------------------------------------------------
cv_knn(x_mat, dise_vec_na, EOC$V, type = "mahala", plot = TRUE)

