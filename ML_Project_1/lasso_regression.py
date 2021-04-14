import sys
#https://www.pluralsight.com/guides/linear-lasso-ridge-regression-scikit-learn
import pandas as pd
import numpy as np
from numpy import arange

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import Lasso
#from sklearn import metrics

encoded_genotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR_bonf_corr_value/PARSED_GENOTYPES_OUTPUT_DIR/encoded_genotypes.tsv"

# The size of the test set for x and y.
test_size = 0.30

# Number of k for k-fold cross validation
k = 10

# Number of repeats.
num_repeats = 3

# Load the the encoded_genotypes_infile to a pandas dataframe encoded_genotypes_df
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
x = encoded_genotypes_df.drop(['genotype_id','phenotype'], axis=1).to_numpy()

# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
y = encoded_genotypes_df['phenotype'].to_numpy(dtype='float32')

print(x)

sys.exit()

print(y)

sys.exit()
# Use the lasso regression (L1 regularization) model.
model = Lasso()

# Define the model cross validation using k = 10 for the evaluation method.
cv = RepeatedKFold(n_splits=k, n_repeats=num_repeats, random_state=1)

# Obtain a grid of alpha values to test to optimize the
grid = dict()
grid['alpha'] = arange(0, 1, 0.01)

# Perform a grid search using the GridSearchCV to fine tune the hyperparameter alpha.
search = GridSearchCV(model, grid, scoring='r2', cv=cv, n_jobs=-1)

# Perform the search and fit the using the model.
results = search.fit(X, y)

# Summarize the model.
print('CV Mean R^2: %.3f' % results.best_score_)
print('Config: %s' % results.best_params_)


# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets. Shuffle contents before splitting.
#x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=float(test_size), shuffle=True)




#model_lasso = Lasso(alpha=0.01)
#model_lasso.fit(X_train, y_train)
#pred_train_lasso= model_lasso.predict(X_train)
#print(np.sqrt(mean_squared_error(y_train,pred_train_lasso)))
#print(r2_score(y_train, pred_train_lasso))
#
#pred_test_lasso= model_lasso.predict(X_test)
#print(np.sqrt(mean_squared_error(y_test,pred_test_lasso)))
#print(r2_score(y_test, pred_test_lasso))
