library(glmnet)
library(MASS)
library(survival)
library(rms)

data <- read.table('CGGA.txt',head=T,row.names=1)
data <- na.omit(data)


yavars<-names(data) %in% c("OS", "Censor")
yavars
CandidateVariables <- as.data.frame(data[!yavars],row.names = NULL)

# Make data frame into matrix
tmp.y <- Surv(data$OS,data$Censor)

tmp.x <- model.matrix(~.,data=CandidateVariables)

# Fit the model
model.lasso <-  glmnet(tmp.x, tmp.y, family="cox", nlambda=50, alpha=1, standardize=TRUE)
plot(model.lasso,xvar="lambda",label=TRUE)

# find the optimal model via cross-validation
cv.model <- cv.glmnet(tmp.x, tmp.y, family="cox", nlambda=50, alpha=1, standardize=TRUE)
plot(cv.model)
cv.model$lambda.min
coef(cv.model, s=cv.model$lambda.min) 

# increase lambda for further shrinkage
cv.model$lambda.1se
Coefficients =coef(cv.model, s=cv.model$lambda.1se) 
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
rowname=row.names(Coefficients)[Active.Index]
result2=cbind(rowname,Active.Coefficients)
write.csv(result2,"LASSO.csv")


