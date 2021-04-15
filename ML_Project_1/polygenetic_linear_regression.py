import sys
#https://www.pluralsight.com/guides/linear-lasso-ridge-regression-scikit-learn
import pandas as pd
import numpy as np
from numpy import arange

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

def mape(y_actual,y_predicted):
    mape = np.mean(np.abs((y_actual - y_predicted)/y_actual))*100
    return mape
    
encoded_genotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR_bonf_corr_value/PARSED_GENOTYPES_OUTPUT_DIR/encoded_genotypes.tsv"

# The size of the train set for x and y.
train_size = 0.70

# The size of the test set for x and y.
test_size = 0.30


# Load the the encoded_genotypes_infile to a pandas dataframe encoded_genotypes_df
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
x = encoded_genotypes_df.drop(['genotype_id','phenotype'], axis=1).to_numpy()

# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
y = encoded_genotypes_df['phenotype'].to_numpy(dtype='float32')


# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets.
x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=float(train_size), test_size=float(test_size))

linear_model = LinearRegression()


model = LinearRegression()
model.fit(x_train, y_train)
pred_train = model.predict(x_train)
mape_train = mape(y_train,pred_train)
print("mape_train value: ", mape_train)
train_accuracy = (100 - mape_train)
print('Train accuracy of Linear Regression: {:0.2f}%.'.format(train_accuracy))
print("mean_squared_error train score: ", mean_squared_error(y_train,pred_train))
print("sqrt mean_squared_error train score: ", np.sqrt(mean_squared_error(y_train,pred_train)))
print("r2_train score: ", r2_score(y_train, pred_train))

pred_test = model.predict(x_test)
mape_test = mape(y_test,pred_test)

print("mape_test value: ", mape_test)
test_accuracy = (100 - mape_test)
print('Test accuracy of Linear Regression: {:0.2f}%.'.format(test_accuracy))
print("mean_squared_error test score: ", mean_squared_error(y_test,pred_test))
print("sqrt mean_squared_error test score: ", np.sqrt(mean_squared_error(y_test,pred_test)))
print("r2_test score: ", r2_score(y_test, pred_test))
