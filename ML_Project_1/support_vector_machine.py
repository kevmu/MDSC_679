import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import metrics

encoded_genotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_output_dir/parsed_genotypes_output_dir/encoded_genotypes.tsv"
phenotypes_infile = "/Users/kevin.muirhead/Desktop/GWAS_output_dir/parsed_genotypes_output_dir/phenotypes.tsv"

test_size = 0.30

### Might need to reformat this file so that we can merge 
encoded_genotypes_df = pd.read_table(encoded_genotypes_infile, delimiter='\t')
phenotypes_df = pd.read_table(phenotypes_infile, delimiter='\t')
print(phenotypes_df)
sys.exit()

# Convert the genotypes x pandas dataframe to numpy array for input into the machine learning algorithm.
#x = dataframe_x.drop('MARKER_ID', axis=1).to_numpy()
merged_genotypes_phenotypes_df=pd.merge(encoded_genotypes_df,phenotypes_df,on=['genotype_id'])

print(merged_genotypes_phenotypes_df)
sys.exit()

# Convert the phenotypes y pandas dataframe to numpy array for input into the machine learning algorithm.
#y = dataframe_y['phenotype'].to_numpy(dtype='float32')
print(y)
#sys.exit()

# Split the X and Y numpy arrays into testing and training. Using training, validation and 30% for test sets. Shuffle contents before splitting.
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=float(test_size), shuffle=True)

sys.exit()
# Use support vector classifier (SVM) with Radial Basis Function (RBF) kernel.
svclassifier = SVC(kernel='rbf', random_state=42)
svclassifier.fit(x_train, y_train)

#https://www.google.com/search?q=train+test+split+sklearn&sxsrf=ALeKk03wUySMMs70aq6eZnyVFZUo4NJZ2A%3A1617769795454&source=hp&ei=QzVtYKiSGZjM-gS3mr6YAg&iflsig=AINFCbYAAAAAYG1DU4GKDf3K0C36p6IZDPLakk2fNT9o&oq=train&gs_lcp=Cgdnd3Mtd2l6EAMYATIECCMQJzIECAAQQzIKCAAQsQMQgwEQQzIECAAQQzIFCAAQsQMyBQguELEDMgUILhCxAzIFCC4QsQMyCAguEMcBEK8BMgUILhCxAzoHCCMQ6gIQJzoICAAQsQMQgwE6CwguELEDEMcBEKMCOgoIABCxAxCDARAKOgIIAFCDSFivTWD5YGgBcAB4AIABaIgBvQOSAQM0LjGYAQCgAQGqAQdnd3Mtd2l6sAEK&sclient=gws-wiz
#svclassifier.score(X_test, y_test)

#scores = cross_val_score(svclassifier, x, y, cv=5, scoring='f1_macro')


#y_pred = svclassifier.predict(x_test)

#print(confusion_matrix(y_test,y_pred))
#print(classification_report(y_test,y_pred))


