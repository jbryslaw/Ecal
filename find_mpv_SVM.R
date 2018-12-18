#find mpv

all_data = read.table(file="out_all.txt",header=TRUE) #all_data = read.table(file="small_test_all.txt",header=TRUE)
attach(all_data)

# SVM is in e1071 library
library(e1071)
#first prepare matrices to pass to SVM:
#0) Matrix with predictors associated with training data:

# m_train = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6, d2ydx2_7)[1:100000,]
# m_test  = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6, d2ydx2_7)[100001:200000,]

# divided by 10
# m_train = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6)[1:10000,]
# m_test  = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6)[10001:20000,]

# m_train = cbind(Min1, Min2, Max0, Max1, Max2, Max3, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4)[1:10000,]     
# m_test  = cbind(Min1, Min2, Max0, Max1, Max2, Max3, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4)[10001:20000,]

m_train = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6)[1:10000,]     
m_test  = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6)[10001:20000,]


# m_train = cbind(Min1, Min2, Max0, Max1, Max2, Max3, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4)[1:100000,]     
# m_test  = cbind(Min1, Min2, Max0, Max1, Max2, Max3, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4)[100001:200000,]

#  divided by 10
v_train_response = MPV[1:10000]     
v_test_response = MPV[10001:20000]

# v_train_response = MPV[1:100000]     
#  v_test_response = MPV[100001:200000]

dat.train = data.frame(x=m_train,y=v_train_response)
dat.test  = data.frame(x=m_test,y=v_test_response)

set.seed(1)

#fit_svm=svm(y~., data=dat.train, kernel="linear", cost=10)
fit_svm=svm(y~., data=dat.train, kernel="radial", cost=1)
#fit_svm=svm(y~., data=dat.train, kernel="polynomial", cost=10)# 0.432 for 1000 evts
#tune to get best cost
#fit_tune=tune(svm,y~., data=dat.train, kernel="polynomial", ranges=list(cost=cbind(0.1 ,1 ,10 ,100 ,1000) ) )
print(summary(fit_svm))
print(fit_svm$fitted)

print("training success rate:")
print(mean( abs(v_train_response-fit_svm$fitted) <=1.0))
hist(abs(v_train_response-fit_svm$fitted))

pred_svm=predict(fit_svm, newdata=dat.test)

print(mean( abs(v_test_response-pred_svm) <=1.0))
pdf("out_1_2_derivs_SVMR1.pdf", useDingbats=FALSE)
hist(abs(v_test_response-pred_svm))
dev.off()
