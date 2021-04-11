import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.model_selection import cross_val_score
from sklearn import metrics

encoded_genotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_output_dir/parsed_genotypes_output_dir/encoded_genotypes.tsv"

test_size = 0.30

# Load the the encoded_genotypes_infile to a pandas dataframe encoded_genotypes_df
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
x = encoded_genotypes_df.drop(['genotype_id','phenotype'], axis=1).to_numpy()


# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
y = encoded_genotypes_df['phenotype'].to_numpy(dtype='float32')

print(x)
print(y)
#sys.exit()






# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets. Shuffle contents before splitting.
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=float(test_size), shuffle=True)

#sys.exit()
# Use support vector machine (SVM) regression model with Radial Basis Function (RBF) kernel.
svr_rbf = SVR(kernel='rbf')
#svr_rbf = SVR(kernel='rbf')
svr_rbf.fit(x_train, y_train)

#https://www.google.com/search?q=train+test+split+sklearn&sxsrf=ALeKk03wUySMMs70aq6eZnyVFZUo4NJZ2A%3A1617769795454&source=hp&ei=QzVtYKiSGZjM-gS3mr6YAg&iflsig=AINFCbYAAAAAYG1DU4GKDf3K0C36p6IZDPLakk2fNT9o&oq=train&gs_lcp=Cgdnd3Mtd2l6EAMYATIECCMQJzIECAAQQzIKCAAQsQMQgwEQQzIECAAQQzIFCAAQsQMyBQguELEDMgUILhCxAzIFCC4QsQMyCAguEMcBEK8BMgUILhCxAzoHCCMQ6gIQJzoICAAQsQMQgwE6CwguELEDEMcBEKMCOgoIABCxAxCDARAKOgIIAFCDSFivTWD5YGgBcAB4AIABaIgBvQOSAQM0LjGYAQCgAQGqAQdnd3Mtd2l6sAEK&sclient=gws-wiz
#print(svr_rbf.score(x_test, y_test))

#sys.exit()
y_pred = svr_rbf.predict(x_test)

print(y_pred)

from sklearn.metrics import mean_absolute_error

print(mean_absolute_error(y_test, y_pred))
sys.exit()
roc_auc_score(y, y_score)

