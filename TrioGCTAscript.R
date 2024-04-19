library(genio) # Read plink format
library(OpenMx) # Fit model
library(gaston) # Compute GRM

# Read data
dat = read.bed.matrix("MCS_topmed")

# Compute GRM
A = GRM(dat)

load("~/trioIDs.RData") # A list with the row numbers that contain mothers, fathers, and offspring
load("~/yX.RData") # A matrix with my phenotype in the first column and five columns of covariates

# Subset GRM blocks
mid = trioIDs[[1]]
pid = trioIDs[[2]]
oid = trioIDs[[3]]
Amm = A[mid, mid]
App = A[pid, pid]
Aoo = A[oid, oid]
Dpm = A[pid, mid] + A[mid, pid]
Dom = A[oid, mid] + A[mid, oid]
Dop = A[oid, pid] + A[pid, oid]


# Set up model

K = 2543 # Number of trios
# Set up the trio model
modmx = mxModel("trio",
                mxMatrix("Lo", 3, 3, T, sqrt(diag(var(yX[, 1]) / 4, 3)),
                         c("l11", "l21", "l31", "l22", "l32", "l33"), name = "L"), # Cholesky factor of genetic covariance matrix
                mxAlgebra(L %*% t(L), name = "Sg"),
                mxMatrix("Fu", 1, 1, T, var(yX[, 1]) / 4, "se", name = "Se"), # Residual variance
                mxMatrix("Sy", K, K, F, Amm, name = "Amm"),
                mxMatrix("Sy", K, K, F, App, name = "App"),
                mxMatrix("Sy", K, K, F, Aoo, name = "Aoo"),
                mxMatrix("Sy", K, K, F, Dpm, name = "Dpm"),
                mxMatrix("Sy", K, K, F, Dom, name = "Dom"),
                mxMatrix("Sy", K, K, F, Dop, name = "Dop"),
                mxMatrix("Id", K, K, name = "I"),
                mxData(yX, "raw", sort = F),
                # Model implied covariance
                mxAlgebra(Sg[1, 1] * Amm + Sg[2, 2] * App + Sg[3, 3] * Aoo +
                            Sg[2, 1] * Dpm + Sg[3, 1] * Dom + Sg[3, 2] * Dop +
                            se * I, name = "V"),
                mxExpectationGREML("V", dataset.is.yX = T),
                mxFitFunctionGREML())

# Fit the model
modmx = mxRun(modmx)

# Look at results
summary(modmx)
mxEval(Sg, modmx) # Genetic (co)variances
mxEval(se, modmx)

#Set up the no parents model
modmx2 = mxModel("noparent",
                mxMatrix("Lo", 1, 1, T, sqrt(var(yX[, 1]) / 2), "l", name = "L"), # Cholesky factor of genetic covariance matrix
                mxAlgebra(L %*% t(L), name = "Sg"),
                mxMatrix("Fu", 1, 1, T, var(yX[, 1]) / 2, "se", name = "Se"), # Residual variance
                mxMatrix("Sy", K, K, F, Aoo, name = "Aoo"),
                mxMatrix("Id", K, K, name = "I"),
                mxData(yX, "raw", sort = F),
                # Model implied covariance
                mxAlgebra( Sg[1, 1] * Aoo +
                            se * I, name = "V"),
                mxExpectationGREML("V", dataset.is.yX = T),
                mxFitFunctionGREML())

# Fit the model
modmx2 = mxRun(modmx2)
summary(modmx2)
mxEval(Sg, modmx2)
mxEval(se, modmx2)


#Set maternal and paternal effects to equal for the combined model
App = (Amm + App + Dpm) / 2
Dop = (Dom + Dop) / 2

#Set up the combined model
modmx3 = mxModel("combined",
                 mxMatrix("Lo", 2, 2, T, sqrt(diag(var(yX[, 1]) / 3, 2)),
                          c("l11", "l21", "l22"), name = "L"), # Cholesky factor of genetic covariance matrix
                 mxAlgebra(L %*% t(L), name = "Sg"),
                 mxMatrix("Fu", 1, 1, T, var(yX[, 1]) / 3, "se", name = "Se"), # Residual variance
                 mxMatrix("Sy", K, K, F, App, name = "App"),
                 mxMatrix("Sy", K, K, F, Aoo, name = "Aoo"),
                 mxMatrix("Sy", K, K, F, Dop, name = "Dop"),
                 mxMatrix("Id", K, K, name = "I"),
                 mxData(yX, "raw", sort = F),
                 # Model implied covariance
                 mxAlgebra(Sg[1, 1] * App + Sg[2, 1] * Dop + Sg[2, 2] * Aoo + se * I, name = "V"),
                 mxExpectationGREML("V", dataset.is.yX = T),
                 mxFitFunctionGREML())

modmx3 = mxRun(modmx3)
summary(modmx3)
mxEval(Sg, modmx3)
mxEval(se, modmx3)

#Set up the null model
modmx4 = mxModel("null",
                 mxMatrix("Fu", 1, 1, T, var(yX[, 1]) / 2, "se", name = "Se"), # Residual variance
                 mxMatrix("Id", K, K, name = "I"),
                 mxData(yX, "raw", sort = F),
                 # Model implied covariance
                 mxAlgebra(se * I, name = "V"),
                 mxExpectationGREML("V", dataset.is.yX = T),
                 mxFitFunctionGREML())

# Fit the model
modmx4 = mxRun(modmx4)
summary(modmx4)
mxEval(Sg, modmx4)
mxEval(se, modmx4)

#Compare Models
mxCompare(modmx,modmx3) #Trio to combined
mxCompare(modmx3,modmx2) #Combined to no parents
mxCompare(modmx2,modmx4) #No parents to null
