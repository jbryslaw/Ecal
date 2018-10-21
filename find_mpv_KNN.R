#find mpv

all_data = read.table(file="out_all.txt",header=TRUE)
attach(all_data)

# KNN is in class library
library(class)
#first prepare matrices to pass to KNN:
#0) Matrices with predictors associated with training data:

#just use min, max and dydx
m_train = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4)[1:100000,]
m_test  = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4)[100001:200000,]

# use min, max, dydx and d2ydx2
# m_train = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6, d2ydx2_7)[1:100000,]	
# m_test  = cbind(Min0, Min1, Min2, Min3, Min4, Max0, Max1, Max2, Max3, Max4, dydx0, dydx1, dydx2, dydx3, dydx4, d2ydx2_0, d2ydx2_1, d2ydx2_2, d2ydx2_3, d2ydx2_4, d2ydx2_5, d2ydx2_6, d2ydx2_7)[100001:200000,]

#1) Fill vectors with response values
v_train_response = MPV[1:100000]     

# Set Random Seed`
set.seed(1)

#2) Perform KNN
fit_knn=knn(m_train,m_test,v_train_response,k=100)

print("results from kNN, k=100:")
# print(fit_knn)
# print("")
# print(v_test_response)
#print(table(fit_knn,v_test_response))

# Calculate how many of the test predictions are within 1 of the true value
print("mean fraction of correct entries (near 1):")
print(mean( abs(v_test_response-as.numeric(paste(fit_knn)))<=1.0))

# Histogram how many of the test predictions are within 1 of the true value
#pdf("out_1_derivs_K100.pdf", useDingbats=FALSE)
hist(abs(as.numeric(paste(fit_knn))-as.numeric(v_test_response)))
#dev.off()
