from itertools import combinations, chain
import csv
import os
import sys

'''

Name: Kevin Muirhead
UCID #: 00502756

MDSC 679: ML Project #2: Implementing the Apriori algorithm

Class Apriori

Performs the AprioriTID algorithm from the following paper;

"An Improved Apriori Algorithm For Association Rules.:, Al-Maolegi, Mohammed & Arkok, Bassam. (2014). International Journal on Natural Language Computing. 3. 10.5121/ijnlc.2014.3103.

Prints out a file that contains the most frequent itemsets and their corresponding association rule metrics; support, confidence, lift, leverage and conviction.

The association rule metric formulas were derived from the following article;

http://rasbt.github.io/mlxtend/user_guide/frequent_patterns/association_rules/

The description for each association rule metric method in this class was adapted from this article.

'''
class Apriori:
    
    '''
    Initialize the Apriori Class.
    
    Input:
    
    database_infile - The database file as input. Contains two columns. The first column contains the transaction id and the second column the transaction entry. The transaction id can consist of any of the following integers, strings or characters. The transaction value entry contains a list of items delimited with a comma ',' character. The items can be any string. Can be integers, floats, strings or characters as long as commas ',' are used to separate items.
    
    association_rule_metrics_outfile -
    
    min_support_count - The minimum support count to consider items for prunning and retaining.
    
    min_confidence - The minimum confidence to consider frequent itemsets
    
    
    
    '''
    def __init__(self, database_infile, association_rule_metrics_outfile, min_support_count, min_confidence):
        
        # The minimum support count to retain items.
        self.min_support_count = int(min_support_count)
        
        # The minimum confidence to retain items.
        self.min_confidence = float(min_confidence)
        
        # The transaction database data structure.
        self.transaction_database = self.load_transaction_database(database_infile)
        
        # Run the apriori algorithm and return the frequent_itemsets data structure.
        self.frequent_itemsets = self.apriori()
        
        # Calculate and print association rules to the association rule metrics output file.
        self.print_association_rules(association_rule_metrics_outfile)

        
    '''
    Load the transaction database using the database infile.
    
    load_transaction_database
    
    Input:
    
    database_infile
        
    Output:
    
    transaction_database - The transaction database data structure. Transaction ID as the key and the Transaction value entry as the value.
    
    '''
    def load_transaction_database(self, database_infile):

        # Counter for header and entry.
        row_counter = 0

        # Transaction database data structure.
        transaction_database = {}
        
        with open(database_infile, 'r') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
            
            # Iterate over each row in the file.
            for row in csv_reader:
            
                #print(row_counter)
                
                # If the line is not the header.
                if(row_counter != 0):
                
                    # Ensure that there are only two items per line.
                    if(len(row) == 2):
                    
                        transaction_id = row[0]
                        transaction_values = row[1]
                        
                        # Remove any spaces between items.
                        transaction_values.replace(" ", "")
                        
                        #print(transaction_id, transaction_values)
#                    else:
#                        print("There are more than one column in this database file")
#                        print(" ".join(["Database input file:", database_infile]))
#                        print(" ".join(["Entry without two columns:", str(row)]))
#                        print(" ".join(["Number of columns:", str(len(row))]))
#                        print(" ".join(["Row number:", str(row_counter)]))
#                        print("Program has been terminated! Please look at the format of the database input file and make sure there are two columns delimited by tabs.")
#                        sys.exit(1)
                    
                    transaction_database[transaction_id] = transaction_values
                    
                row_counter += 1
                
        return transaction_database
                


    '''
    
    Method generate_candidate_set1()

    Output:

    large_1_itemset - Returns the large 1-Itemset of the transaction database.
        
    '''
    def generate_candidate_set1(self):
        
        # The large_1_itemset support dictionary data structure.
        large_1_support_count = {}
        
        # The large_1_itemset data list structure.
        large_1_itemset = []
        
        # The large_1_transaction_ids dictionary data structure.
        large_1_transaction_ids = {}
        
        # Iterate over each transaction in the database.
        for transaction_id in self.transaction_database:
            #print(transaction_id)
            #sys.exit()
            #print(transaction_database[transaction])
            items = self.transaction_database[transaction_id].split(",")
            
            # Iterate over the items in the transaction.
            for item in items:
                #print(item)
                
                # If item does not exist in L1 then set the support counter to 1.
                # Otherwise increment the item by 1.
                if(not(item in large_1_support_count)):
                    large_1_support_count[item] = 1
                    
                elif(item in large_1_support_count):
                    large_1_support_count[item] += 1
                    
                if(not(item in large_1_transaction_ids)):
                    large_1_transaction_ids[item] = []
                    large_1_transaction_ids[item].append(transaction_id)
                elif(item in large_1_transaction_ids):
                    large_1_transaction_ids[item].append(transaction_id)
                            
        
        # Get the keys of the large_1_support_count data structure and make the sorted large_1_itemset
        large_1_itemset = list(sorted(large_1_support_count.keys()))
                
        # The item_list data structure for comparing items using support counts.
        item_list = []
        
        # For each item in large_1_itemset compare the support_count to the min_support_count.
        for item in large_1_itemset:
            #print(item)
            
            # The support_count of the item in large_1_itemset.
            support_count = large_1_support_count[item]
            
            #print(item)
            #print(support_count)
            
            # If item has a support_count less than the min_support_count then append item to item_list.
            if(support_count < self.min_support_count):
                #print("support_count < self.min_support_count")
                #print(str(support_count) + " < " + str(self.min_support_count))
                item_list.append(item)
        
        #print(item_list)
        
        # Remove items from large_1_itemset that are in the item_list.
        if(len(item_list) > 1):
            for item in sorted(item_list):
                large_1_itemset.remove(item)
            
        del item_list
            
        return (large_1_itemset, large_1_support_count, large_1_transaction_ids)


    '''
    apriori_gen - Apriori candidate generation step.  Join step

    '''
    def apriori_gen(self, prev_candidate_itemset, k):

        # The one itemset of the prev_candidate_itemset (previous candidate itemset) for making combinations for candidates.
        previous_1_itemset = []
        
        # If prev_candidate_itemset is the large_1_itemset (Large 1-Itemset).
        if((k - 1) == 1):
        
            previous_1_itemset = prev_candidate_itemset
            
        # Else prev_candidate_itemset is the Lk_itemset (Large k-Itemset).
        elif((k - 1) > 1):
            
            #print(prev_candidate_itemset)
            
            # A dictionary datastructure to get the set of previous_1_itemset
            previous_1_itemset_dict = {}
            
            # Iterate through all the items in the prev_candidate_itemset.
            for prev_k_item in prev_candidate_itemset:
                prev_k_item = list(prev_k_item)
                #print(prev_k_item)
                
                # Iterate through each item in the prev_k_item list.
                for item in prev_k_item:
                    previous_1_itemset_dict[item] = item
            
    #        print(sorted(list(previous_1_itemset_dict.keys())))
            previous_1_itemset = sorted(list(previous_1_itemset_dict.keys()))
    #        print(previous_1_itemset)
            #sys.exit()
        
        # Generate all possible combinations from the previous_1_itemset.
        candidate_k_itemset = combinations(previous_1_itemset, k)
        #print(list(candidate_k_itemset))
        
        # Remove k-items from candidate_k_itemset so that we have the new candidate_k_itemset of items.
        new_candidate_k_itemset = []
        for item_set in list(candidate_k_itemset):
            #print(item_set)
            num_equal_items = 0
            for p_item in item_set:
                #print(p_item)
                #sys.exit()
                for q_item in previous_1_itemset:
                    #print(q_item)
                    #sys.exit()
                    if(p_item == q_item):
                        #print(p_item,q_item)
                        
                        num_equal_items = num_equal_items + 1
                
                # If the number of items in p_item equal to q_item is equal to k.
                if(num_equal_items == k):
                    #print(item_set)
                    new_candidate_k_itemset.append(item_set)
                    break
                    
        #print(new_candidate_k_itemset)
        del candidate_k_itemset

        return(new_candidate_k_itemset)


    def prune_k_itemset(self, candidate_support_count, candidate_transaction_ids):
        #print(candidate_support_count)
        
        # The large-k itemset data structure.
        large_k_itemset = []
        
        # The large-k support count data structure.
        large_k_support_count = {}
            
        # The large-k transaction ids data structure.
        large_k_transaction_ids = {}
        
        # Iterate over all the candidate counts.
        for k_item_key in candidate_support_count:
        
            # The support_count of the item in large_1_itemset.
            support_count = candidate_support_count[k_item_key]
            
            #print(k_item_key)
            #print(support_count)
            
            #sys.exit()
            
            # If item has a support_count less than the min_support_count then append item to item_list.
            if(support_count >= self.min_support_count):
                
                #print(k_item_key)
                #sys.exit()
                k_item = tuple(k_item_key.split(","))
                #sys.exit()
                large_k_itemset.append(k_item)
                large_k_support_count[k_item_key] = candidate_support_count[k_item_key]
                large_k_transaction_ids[k_item_key] = candidate_transaction_ids[k_item_key]
           
        # Delete the candidate_support_count data structure.
        del candidate_support_count
        del candidate_transaction_ids
        
        return (large_k_itemset, large_k_support_count, large_k_transaction_ids)
        
        
    '''
    Function apriori

    Description:

    The Apriori algorithm

    Output:

    frequent_itemsets - The frequent itemset dictionary data structure containing the list of frequent itemsets and metadata for each value of k. The frequent_itemsets data structure contains a dictionary of frequent itemsets, a dictionary containing frequent itemset support counts, and a dictionary of itemsets containing a list of transaction ids that contain the itemset.


    '''
    def apriori(self):

        # The frequent_itemsets dictionary data structure.
        frequent_itemsets = {}
        
        # Generate candidate set 1 and get the large 1-itemset and the large 1 itemset support counts as well as the transaction ids for each item in the large 1-itemset.
        (large_1_itemset, large_1_support_count, large_1_transaction_ids) = self.generate_candidate_set1()

        # Set k = 1 for 1-Itemset candidate generation.
        k = 1
        
        # Generate frequent_itemsets at k = 1. Add the data structure for large_1_itemset, large_1_support_count, large_1_transaction_ids
        itemset_dict = {}
        itemset_dict['itemset_data'] = large_1_itemset
        itemset_dict['itemset_support_count'] = large_1_support_count
        itemset_dict['itemset_transaction_ids'] = large_1_transaction_ids
        frequent_itemsets[str(k)] = itemset_dict
        
        # Set the current k_itemset to the large_1_itemset. Going forward this will be the previous k_itemset.
        previous_k_itemset = large_1_itemset
        
        # Set k = 2 to start the iteration step.
        k = 2
        
        # Iterate over the Large k-itemsets until the previous_k_itemset is empty. If previous_k_itemset is empty then terminate the while loop.
        while(len(previous_k_itemset) > 0):
            
            # Generate new apriori candidates.
            candidate_itemset_k = self.apriori_gen(previous_k_itemset, k)

            # The candidate in transaction counter.
            candidate_support_count = {}
            
            # The candidate in transaction ids dictionary data structure.
            candidate_transaction_ids = {}

            # Set all the k_items in candidate_itemset_k to zero in the candidate_support_count.
            for k_item in candidate_itemset_k:
                k_item = list(k_item)
                k_item_key = ",".join(k_item)
                #print(k_item_key)
                candidate_support_count[k_item_key] = 0

            #print(candidate_support_count)
            #sys.exit()

            # Iterate over each transaction in the database.
            for transaction_id in self.transaction_database:
                #print(transaction_id)
                #print(transaction_database[transaction_id])
                
                # Split the transaction into a list of items.
                k_itemA = self.transaction_database[transaction_id].split(",")
                
                # Convert the k_itemA list to a set. So we can figure out if there are subsets in candidate_itemset_k
                k_itemA = set(k_itemA)
                for k_itemB in candidate_itemset_k:
                    #print(k_itemB)
                    k_itemB = set(k_itemB)
                    #print(k_itemB)
                    # If k_itemB is a subset of k_itemA then increment
                    if(k_itemB.issubset(k_itemA)):
                        #print(k_itemA, k_itemB)
                        #print(k_itemB.issubset(k_itemA))
                        k_itemB = sorted(list(k_itemB))
                        
                        # Join the sorted k_itemB list into a key for the candidate_support_count to increment the counter.
                        k_itemB_key = ",".join(k_itemB)
                        #print(k_itemB_key)
                        
                        # Increment candidate_support_count at k_itemB_key because k_itemB_key is in the transaction.
                        candidate_support_count[k_itemB_key] = candidate_support_count[k_itemB_key] + 1
                        if(not(k_itemB_key in candidate_transaction_ids)):
                            candidate_transaction_ids[k_itemB_key] = []
                            candidate_transaction_ids[k_itemB_key].append(transaction_id)
                        elif(k_itemB_key in candidate_transaction_ids):
                            candidate_transaction_ids[k_itemB_key].append(transaction_id)
                    
            #print(candidate_support_count)

            #print(candidate_transaction_ids)
            #sys.exit()
            
            # Prune items from candidate k-itemset that does not meet the following threshold support_count >= min_support_count and return the large_k_itemset, large_k_support_count and large_k_transaction_ids.
            (large_k_itemset, large_k_support_count, large_k_transaction_ids) = self.prune_k_itemset(candidate_support_count, candidate_transaction_ids)

            # Generate frequent_itemsets at k if large_k_itemset is not empty (length is equal to zero).
            if(len(large_k_itemset) > 0):
                itemset_dict = {}
                itemset_dict['itemset_data'] = large_k_itemset
                itemset_dict['itemset_support_count'] = large_k_support_count
                itemset_dict['itemset_transaction_ids'] = large_k_transaction_ids
                frequent_itemsets[str(k)] = itemset_dict
        
            k = k + 1
            
            # Set the previous_k_itemset to the large_k_itemset for next iteration.
            previous_k_itemset = large_k_itemset
            
            #print(frequent_itemsets)
            
        return(frequent_itemsets)

    '''
    
    '''
    def support(self, item_support_count, num_transactions):
        support = (float(item_support_count)/float(num_transactions))
        return support
    
    '''
    
    '''
    def confidence(self, itemA_support, itemB_support):
        
        # confidence = (support(itemA U itemB))/(support(itemB))
        
        # Calculation example from https://www.geeksforgeeks.org/apriori-algorithm/
        # [I1^I2]=>[I3] //confidence = sup(I1^I2^I3)/sup(I1^I2) = 2/4*100=50%
        # [I1^I3]=>[I2] //confidence = sup(I1^I2^I3)/sup(I1^I3) = 2/4*100=50%
        # [I2^I3]=>[I1] //confidence = sup(I1^I2^I3)/sup(I2^I3) = 2/4*100=50%
        # [I1]=>[I2^I3] //confidence = sup(I1^I2^I3)/sup(I1) = 2/6*100=33%
        # [I2]=>[I1^I3] //confidence = sup(I1^I2^I3)/sup(I2) = 2/7*100=28%
        # [I3]=>[I1^I2] //confidence = sup(I1^I2^I3)/sup(I3) = 2/6*100=33%
        
        # Initialize confidence variable to zero float to start.
        confidence = 0.00
        
        # Check to make sure that we aren't dividing by zero.
        if(float(itemB_support) != 0):
            confidence = float(itemA_support)/float(itemB_support)
        else: # If zero then throw error message and terminate.
            print("Cannot divide by zero!")
            print(" ".join(["itemB_support", "=", itemB_support]))
            sys.exit(1)
        
        return confidence
        
    '''
    
    '''
    def lift(self, confidence, itemB_support):
        # lift = confidence(itemA U itemB)/support(itemB)
        lift = 0.00

        # Check to make sure that we aren't dividing by zero.
        if(float(itemB_support) != 0):
            lift = float(confidence)/float(itemB_support)
        else: # If zero then set lift to positive infinity.
            lift = float('inf')
            
        return lift
       
    '''

    '''
    def leverage(self, itemA_support, itemB_support, subset_support):
    
        # leverage = support(itemA U itemB)/(support(itemA) * support(itemB))
        leverage = float(itemA_support)/(float(subset_support) * float(itemB_support))
        
        return leverage
        
    '''

    '''
    def conviction(self, confidence, itemB_support):
    
        # conviction = (1 - support(itemB))/(1 - confidence(itemA U itemB))
        
        # Initialize the conviction to zero float value.
        conviction = 0.00
        
        # Check to make sure that we aren't dividing by zero.
        if(float(1 - float(confidence)) != 0):
            conviction = float(1 - float(itemB_support))/float(1 - float(confidence))
        else: # If zero then set conviction to positive infinity.
            conviction = float('inf')
        
        return conviction
        
    '''
    
    '''
    def print_association_rules(self, association_rule_metrics_outfile):

        # Writing the association rule metrics to a TSV file for summary of association rules for each item in each k-itemset.
        tsvfile = open(association_rule_metrics_outfile, 'w')
        tsvwriter = csv.writer(tsvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
    
        # Write the header for the association rule metrics TSV file.
        tsvwriter.writerow(['itemA', 'transaction_ids', 'itemA_support', 'itemB', 'itemB_support', 'subset', 'subset_support', 'confidence', 'confidence_rules', 'lift', 'leverage', 'conviction'])
        
        # The association rule metrics data structure to contain all
        association_rule_metrics = []
        
        # Iterate over all the frequent itemsets so that we can calculate the association rule metrics.
        for k in sorted(self.frequent_itemsets, reverse=True):
        
            # The frequent k-itemsets containing itemA.
            frequent_k_itemsets_A = self.frequent_itemsets[k]
            
            # The itemset support count variable containing support counts for itemA.
            itemset_support_count_A = frequent_k_itemsets_A["itemset_support_count"]
            
            # The transaction ids that are associated with this itemset.
            itemset_transaction_ids = frequent_k_itemsets_A["itemset_transaction_ids"]
  
            # Iterate over the itemset_support_count data structure for support_counts for each item in the itemset.
            for itemA in sorted(itemset_support_count_A):
                #print(itemA)
                
                # Get the support count of the item in the itemset.
                itemA_support_count = itemset_support_count_A[itemA]
                #print(support_count)
                
                # Split itemA into individual components so that we can find the powerset of the item.
                itemA_list = itemA.split(",")
                
                # Generate powerset list from the itemset_list so that we can calculate all the assocation rules.
                powerset = list(chain.from_iterable(combinations(itemA_list, r) for r in range(len(itemA_list)+1)))
                
                #print(powerset)
                
                # The association rule metrics data structure.
                association_rule_metrics = {}
                for set in powerset:
                
                    # We only want powerset sets that are of length 1 to less than the size of itemset
                    if(len(set) > 0 and len(set) < len(itemA_list)):
                        set_list_itemB = []
                        for item in set:
                            #print(item)
                            set_list_itemB.append(item)
                        #print(set_list_itemB)
                        
                        # Join the list into a string representation of itemB.
                        itemB = ",".join(set_list_itemB)
                        
                        # The frequent k-itemsets containing itemB.
                        frequent_k_itemsets_B = self.frequent_itemsets[str(len(set_list_itemB))]
                        
                        # The itemset support count variable containing support counts for itemB.
                        itemset_support_count_B = frequent_k_itemsets_B["itemset_support_count"]
                        
                        # Get the support count for itemB.
                        itemB_support_count = itemset_support_count_B[str(itemB)]
                        
                        #print(itemB_support_count)
                        
                        # The number of transactions in the database.
                        num_transactions = len(self.transaction_database.keys())
                        
                        # Calculate the support of itemA.
                        # itemA_support = support(itemA)/number of transactions in the database)
                        itemA_support = self.support(itemA_support_count,num_transactions)
                        
                        # Calculate the support of itemB.
                        # itemB_support = support(itemB)/number of transactions in the database)
                        itemB_support = self.support(itemB_support_count,num_transactions)
                        
                        # confidence = support(itemA U itemB)/support(itemB)
                        confidence = self.confidence(itemA_support, itemB_support)
                        
                        # Getting difference of itemA_list and set_list_itemB to obtain a subset of the two for print out of confidence caluclations.
                        subset_list = [i for i in itemA_list + set_list_itemB if i not in itemA_list or i not in set_list_itemB]

                        # Join the list into a string representation of the subset.
                        subset = ",".join(subset_list)
                        #print(subset)
                        
                        # lift = confidence(itemA U itemB)/support(itemB)
                        lift = self.lift(confidence, itemB_support)
                        
                        # The frequent k-itemsets containing the subset of itemA U itemB.
                        frequent_k_itemsets_subset = self.frequent_itemsets[str(len(subset_list))]

                        # The itemset support count variable containing support counts for the subset of itemA U itemB.
                        itemset_support_count_subset = frequent_k_itemsets_subset["itemset_support_count"]
                        
                        # Get the support count for the subset of itemA U itemB.
                        subset_support_count = itemset_support_count_subset[str(subset)]
                        
                        # Calculate the support of the subset item.
                        subset_support = float(subset_support_count)/float(num_transactions)
                        
                        # leverage = (support(itemA U itemB))/(support(itemA) * support(itemB))
                        leverage = self.leverage(itemA_support, itemB_support, subset_support)
                        
                        # conviction = 1 - support(itemB)/1 - confidence(itemA U itemB)
                        conviction = self.conviction(confidence, itemB_support)
                        
                        if(float(confidence) >= float(self.min_confidence)):
                            confidence_rules = "({itemB} => {subset}), confidence = support({itemA})/support({itemB}) = ({itemA_support}/{itemB_support}) * 100 = {confidence}%".format(itemA=itemA, itemB=itemB, subset=subset, itemA_support="{:.3f}".format(itemA_support),itemB_support="{:.3f}".format(itemB_support), confidence="{:.3f}".format(confidence*100))
                            
                            # Get the list of transaction ids for itemA.
                            transaction_ids = ",".join(sorted(itemset_transaction_ids[itemA]))
                            #print(transaction_ids)
                            #print(itemA)
                            
                            # Print the association rules for each item in each k-itemset.
                            tsvwriter.writerow([itemA, transaction_ids, "{:.3f}".format(itemA_support),itemB, "{:.3f}".format(itemB_support),subset, "{:.3f}".format(subset_support), "{:.3f}".format(confidence),confidence_rules, "{:.3f}".format(lift), "{:.3f}".format(leverage), "{:.3f}".format(conviction)])
                            
        # Close the association rule metrics TSV output file.
        tsvfile.close()


#### MAIN FUNCTION FOR TESTING ####

#output_dir = "/Users/kevin.muirhead/Desktop/macbook_air/MDSC_679/ML_Project_2/test_datasets/association_rules"
#
## Create the output directory if it does not exist.
#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)
#    
## The transaction database input file.
#database_infile = "/Users/kevin.muirhead/Desktop/macbook_air/MDSC_679/ML_Project_2/test_datasets/transaction_database1.tsv"
#
## Get the basename of the transaction database input file.
#basename = os.path.basename(database_infile)
#
## Get the filename without the extension so we can make the association rule metrics output filename.
#filename = os.path.splitext(basename)[0]
#
## The association rule metrics output file.
#association_rule_metrics_outfile = os.path.join(output_dir, "_".join([filename, "association_rule_metrics.tsv"]))
#
## Minimum support count is 2.
#min_support_count = 2
#
## Minimum confidence is 60 %.
#min_confidence = 0.60
#
#Apriori(database_infile, association_rule_metrics_outfile, min_support_count, min_confidence)




