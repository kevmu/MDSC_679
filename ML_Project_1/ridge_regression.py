'''
Name: Kevin Muirhead
UCID#: 00502756

ridge_regression.py - Executes RIDGE Regression (L2 Regularization) on a encoded genotypes file as input.

'''

import os
import sys
import csv
import pandas as pd
import numpy as np
from numpy import arange

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import Ridge

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from math import sqrt

import argparse

parser = argparse.ArgumentParser()

encoded_genotypes_infile = None
output_dir = None

# Usage example
# python ridge_regression.py --encoded_genotypes_infile /Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR_bonf_corr_value/PARSED_GENOTYPES_OUTPUT_DIR/encoded_genotypes.tsv --output_dir /Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR_bonf_corr_value/PARSED_GENOTYPES_OUTPUT_DIR/RIDGE_REGRESSION_METRICS_OUTPUT_DIR

parser = argparse.ArgumentParser(description='Perform RIDGE Regression (L2 Regularization) on a encoded genotypes file as input.' )

parser.add_argument('--encoded_genotypes_infile', action='store', dest='encoded_genotypes_infile',
                    help='The encoded genotypes file as input. (i.e. $HOME/filename.tsv)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='The output directory as input. (i.e. $HOME/output_dir)')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

encoded_genotypes_infile = results.encoded_genotypes_infile
output_dir = results.output_dir

if(encoded_genotypes_infile == None):
    print('\n')
    print('error: please use the --encoded_genotypes_infile option to specify the encoded genotypes file as input. (i.e. filename.tsv)')
    print('encoded_genotypes_infile =' + ' ' + str(encoded_genotypes_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

# Create the output_dir directory if it does not already exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get the basename of the encoded genotypes input file.
basename = os.path.basename(encoded_genotypes_infile)

# Get the filename without the extension so we can make the output filename.
filename = os.path.splitext(basename)[0]

# The Mean Absolute Percentage Error (MAPE) score calculation.
def mape(y_actual,y_predicted):
    mape = np.mean(np.abs((y_actual - y_predicted)/y_actual))*100
    return mape

# The size of the train set for x and y.
train_size = 0.70

# The size of the test set for x and y.
test_size = 0.30

# Number of k for k-fold cross validation
k = 5

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
    
# Find optimal alpha with a grid search.
alphas = arange(0, 1, 0.0001)

print(len(alphas))

# Perform the randomized search to optimize hypertuning parameters.
param_grid = {'alpha': alphas}
grid = GridSearchCV(estimator=model, param_grid=dict(alpha=alphas))
grid.fit(x_train, y_train)

# Summarize the results of the grid search
print(grid.best_score_)
print(grid.best_estimator_.alpha)

best_alpha = grid.best_estimator_.alpha


print("Best alpha: ", best_alpha)

# The metrics output file.
metrics_outfile = os.path.join(output_dir, "_".join([filename, "ridge_regression_metrics.tsv"]))

# Writing the model evaluation metrics to a TSV file for summary of metrics
tsvfile = open(metrics_outfile, 'w')
tsvwriter = csv.writer(tsvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

# Write the header for the model evaluation metrics TSV file.
tsvwriter.writerow(["mean_squared_error_train","mean_squared_error_test","root_mean_squared_error_train","root_mean_squared_error_test","r2_score_train","r2_score_test","mape_train","mape_test","mape_train_accuracy","mape_test_accuracy"])
        
        
# The Ridge Regression (L2 Regularization) model.
model = Ridge(alpha=best_alpha, fit_intercept=True, normalize=False, max_iter=1000)

# Fit the model.
model.fit(x_train, y_train)

# Predict the model.
pred_train = model.predict(x_train)

# Calculate the Mean Absolute Percentage Error (MAPE) train set score.
mape_train = mape(y_train,pred_train)
print("mape_train value: ", mape_train)

# The Mean Absolute Percentage Error (MAPE) train set accuracy score.
mape_train_accuracy = (100 - mape_train)
print("Train accuracy of Ridge Regression: {:0.2f}%.".format(mape_train_accuracy))

# The Mean Squared Error (MSE) train set score.
mean_squared_error_train = mean_squared_error(y_train,pred_train)
print("mean_squared_error train score: ", mean_squared_error_train)

# The Root Mean Squared Error (RMSE) train set score.
root_mean_squared_error_train = np.sqrt(mean_squared_error(y_train,pred_train))
print("sqrt mean_squared_error train score: ", root_mean_squared_error_train)

# The Coefficient of determination R-squared (R2) train set score.
r2_score_train = r2_score(y_train, pred_train)
print("r2_train score: ", r2_score_train)

# Predict the model
pred_test = model.predict(x_test)

# Calculate the Mean Absolute Percentage Error (MAPE) test set score.
mape_test = mape(y_test,pred_test)
print("mape_test value: ", mape_test)

# The Mean Absolute Percentage Error (MAPE) test set score.
mape_test_accuracy = (100 - mape_test)

# The Mean Absolute Percentage Error (MAPE) test set accuracy score.
print("Test accuracy of Ridge Regression: {:0.2f}%.".format(mape_test_accuracy))

# The Mean Squared Error (MSE) test set score.
mean_squared_error_test = mean_squared_error(y_test,pred_test)
print("mean_squared_error test score: ", mean_squared_error_test)

# The Root Mean Squared Error (RMSE) test set score.
root_mean_squared_error_test = np.sqrt(mean_squared_error(y_test,pred_test))
print("sqrt mean_squared_error test score: ", root_mean_squared_error_test)

# The Coefficient of determination R-squared (R2) test set score.
r2_score_test = r2_score(y_test, pred_test)
print("r2_test score: ", r2_score_test)

# Write the evaluation metrics to a file.
tsvwriter.writerow([mean_squared_error_train,mean_squared_error_test,root_mean_squared_error_train,root_mean_squared_error_test,r2_score_train,r2_score_test,mape_train,mape_test,mape_train_accuracy,mape_test_accuracy])
    
    
