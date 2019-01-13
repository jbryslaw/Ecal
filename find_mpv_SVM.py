import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from scipy import stats as scipystats
from sklearn import svm

df_total = pd.read_pickle("all_features.plk")
print(df_total)

# use 1st half of the data as training
#    and 2nd half as test data
i_halfpt = np.floor(len(df_total)/2.).astype(int)

print(df_total.shape)
print(i_halfpt)

df_train = df_total.iloc[0:i_halfpt]
df_test  = df_total.iloc[i_halfpt:len(df_total)]
x_train  = df_total.iloc[0:i_halfpt,0:22]
x_test   = df_total.iloc[i_halfpt:len(df_total),0:22]
y_train  = df_total.iloc[0:i_halfpt,22]
y_test   = df_total.iloc[i_halfpt:len(df_total),22]

# print(df_train)
# print(x_train)

fit_svr = svm.SVR()
fit_svr.fit(x_train,y_train)
y_reco  = fit_svr.predict(x_test)


# fit_knn = libknn.KNeighborsRegressor(n_neighbors=1)
# fit_knn.fit(x_train,y_train)
# y_reco  = (fit_knn.predict(x_test))

print(y_reco)
# print(y_test)
# print(abs(y_reco-y_test))

out_fig = plt.figure(figsize=(16,10),dpi=80)
plt.hist(abs(y_reco-y_test))
print("test sucess rate:")
np.mean(abs(y_reco-y_test)<1.0)

