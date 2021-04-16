import sys
#https://www.pluralsight.com/guides/linear-lasso-ridge-regression-scikit-learn
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

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


# Load the the encoded_genotypes_infile to a pandas dataframe encoded_genotypes_df
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
x = encoded_genotypes_df.drop(['genotype_id','phenotype'], axis=1).to_numpy()

# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
y = encoded_genotypes_df['phenotype'].to_numpy(dtype='float32')


# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets.
x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=float(train_size), test_size=float(test_size))

# The Linear Regression model.
model = LinearRegression()

# Fit the model.
model.fit(x_train, y_train)

# Predict the model.
pred_train = model.predict(x_train)

# Calculate the Mean Absolute Percentage Error (MAPE) train set score.
mape_train = mape(y_train,pred_train)

# The Mean Absolute Percentage Error (MAPE) train set score.
print("mape_train value: ", mape_train)

# The Mean Absolute Percentage Error (MAPE) train set accuracy score.
train_accuracy = (100 - mape_train)
print('Train accuracy of Linear Regression: {:0.2f}%.'.format(train_accuracy))

# The Mean Squared Error (MSE) train set score.
print("mean_squared_error train score: ", mean_squared_error(y_train,pred_train))

# The Root Mean Squared Error (RMSE) train set score.
print("sqrt mean_squared_error train score: ", np.sqrt(mean_squared_error(y_train,pred_train)))

# The Coefficient of determination R-squared (R2) train set score.
print("r2_train score: ", r2_score(y_train, pred_train))

# Predict the model.
pred_test = model.predict(x_test)

# Calculate the Mean Absolute Percentage Error (MAPE) test set score.
mape_test = mape(y_test,pred_test)

# The Mean Absolute Percentage Error (MAPE) test set score.
print("mape_test value: ", mape_test)

# The Mean Absolute Percentage Error (MAPE) test set accuracy score.
test_accuracy = (100 - mape_test)
print('Test accuracy of Linear Regression: {:0.2f}%.'.format(test_accuracy))

# The Mean Squared Error (MSE) test set score.
print("mean_squared_error test score: ", mean_squared_error(y_test,pred_test))

# The Root Mean Squared Error (RMSE) test set score.
print("sqrt mean_squared_error test score: ", np.sqrt(mean_squared_error(y_test,pred_test)))

# The Coefficient of determination R-squared (R2) test set score.
print("r2_test score: ", r2_score(y_test, pred_test))
