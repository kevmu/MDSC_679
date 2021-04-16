import sys
import pandas as pd
import numpy as np
from numpy import arange

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import Ridge

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from math import sqrt

# The Mean Absolute Percentage Error (MAPE) score calculation.
def mape(y_actual,y_predicted):
    mape = np.mean(np.abs((y_actual - y_predicted)/y_actual))*100
    return mape

# The encoded genotypes and filtered phenotypes input file.
encoded_genotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR_bonf_corr_value/PARSED_GENOTYPES_OUTPUT_DIR/encoded_genotypes.tsv"

# The size of the train set for x and y.
train_size = 0.70

# The size of the test set for x and y.
test_size = 0.30

# Number of k for k-fold cross validation
k = 10

# Number of repeats. (Number of iterations).
num_repeats = 10

# Load the the encoded_genotypes_infile to a pandas dataframe encoded_genotypes_df
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
x = encoded_genotypes_df.drop(['genotype_id','phenotype'], axis=1).to_numpy()

# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
y = encoded_genotypes_df['phenotype'].to_numpy(dtype='float32')


# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets.
x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=float(train_size), test_size=float(test_size))

# Train model with default alpha=1
model = Ridge(alpha=1,tol=1e-4, fit_intercept=True, normalize=False, max_iter=1000).fit(x_train, y_train)

# Define the model cross validation using k = 10 for the evaluation method.
cv = RepeatedKFold(n_splits=k, n_repeats=num_repeats, random_state=1)

# Get cross validation scores.
scores = cross_val_score(model,
    x_train,
    y_train,
    cv=cv,
    scoring='r2'
)
    
print('CV Mean: ', np.mean(scores))
print('STD: ', np.std(scores))
    
# find optimal alpha with a randomized search
alpha = arange(0, 8, 0.0001)

print(len(alpha))


param_grid = {'alpha': alpha}
rand_search = RandomizedSearchCV(estimator=model,
                                param_distributions=param_grid,
                               n_iter=10000)

rand_search.fit(x_train, y_train)
best_alpha = rand_search.best_estimator_.alpha


print(rand_search.best_estimator_.alpha)
print(rand_search.best_score_)
print(best_alpha)

# The Ridge Regression (L2 Regularization) model.
model = Ridge(alpha=best_alpha, fit_intercept=True, normalize=False, max_iter=1000)

# Fit the model.
model.fit(x_train, y_train)

# Train the model.
pred_train = model.predict(x_train)

# Calculate the Mean Absolute Percentage Error (MAPE) train set score.
mape_train = mape(y_train,pred_train)
print("mape_train value: ", mape_train)

# The Mean Absolute Percentage Error (MAPE) train set accuracy score.
train_accuracy = (100 - mape_train)
print('Train accuracy of Ridge Regression: {:0.2f}%.'.format(train_accuracy))

# The Mean Squared Error (MSE) train set score.
print("mean_squared_error train score: ", mean_squared_error(y_train,pred_train))

# The Root Mean Squared Error (RMSE) train set score.
print("sqrt mean_squared_error train score: ", np.sqrt(mean_squared_error(y_train,pred_train)))

# The Coefficient of determination R-squared (R2) train set score.
print("r2_train score: ", r2_score(y_train, pred_train))

# Predict the model
pred_test = model.predict(x_test)

# Calculate the Mean Absolute Percentage Error (MAPE) test set score.
mape_test = mape(y_test,pred_test)
print("mape_test value: ", mape_test)

# The Mean Absolute Percentage Error (MAPE) test set score.
test_accuracy = (100 - mape_test)

# The Mean Absolute Percentage Error (MAPE) test set accuracy score.
print('Test accuracy of Ridge Regression: {:0.2f}%.'.format(test_accuracy))

# The Mean Squared Error (MSE) test set score.
print("mean_squared_error test score: ", mean_squared_error(y_test,pred_test))

# The Root Mean Squared Error (RMSE) test set score.
print("sqrt mean_squared_error test score: ", np.sqrt(mean_squared_error(y_test,pred_test)))

# The Coefficient of determination R-squared (R2) test set score.
print("r2_test score: ", r2_score(y_test, pred_test))

