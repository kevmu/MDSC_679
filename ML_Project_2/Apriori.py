from itertools import combinations
import csv
import sys

def load_transaction_database(database_infile):

    # Counter for header and entry.
    row_counter = 0

    # Transaction database data structure.
    transaction_database = {}
    
    with open(database_infile, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        
        # Iterate over each row in the file.
        for row in csv_reader:
        
            print(row_counter)
            
            # If the line is not the header.
            if(row_counter != 0):
            
                # Ensure that there are only two items per line.
                if(len(row) == 2):
                
                    transaction_id = row[0]
                    transaction_values = row[1]
                    
                    # Ensure that transaction_values are sorted.
                    ########
                    
                    # Remove any spaces between items
                    transaction_values.replace(" ", "")
                    
                    #print(transaction_id, transaction_values)
                transaction_database[transaction_id] = transaction_values
                
            row_counter += 1
    return transaction_database
            
#def support(item, num_transactions)
#
#    support = (item.frequency / num_transactions)
#    return support

#def confidence()
#
#    confidence = (support(itemX U itemY))/(support(itemX))
#
#    return condfidence
#
#def lift(itemX, itemY)
#
#    lift = (support(itemX U itemY))/(support(itemX)  * support(itemY))
#    return lift



'''
Function generate_candidate_set1(transaction_database, min_support_count)

Input:

transaction_database - The transaction database data structure.

min_support_count - The minimum support count for excluding items in the 1-Itemset.

Output:

Returns:

L1 - The 1-Itemset of the transaction database.
    
'''
def generate_candidate_set1(transaction_database, min_support_count):
    
    # The large_1_itemset support dictionary data structure.
    large_1_support_count = {}
    
    # The large_1_itemset data list structure.
    large_1_itemset = []
    
    # The large_1_transaction_ids dictionary data structure.
    large_1_transaction_ids = {}
    
    # Iterate over each transaction in the database.
    for transaction_id in transaction_database:
        print(transaction_id)
        #sys.exit()
        #print(transaction_database[transaction])
        items = transaction_database[transaction_id].split(",")
        
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
        if(support_count < min_support_count):
            #print("support_count < min_support_count")
            #print(str(support_count) + " < " + str(min_support_count))
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
def apriori_gen(prev_candidate_itemset, k):

    # The one itemset of the prev_candidate_itemset (previous candidate itemset) for making combinations for candidates.
    previous_1_itemset = []
    
    # If prev_candidate_itemset is the large_1_itemset (Large 1-Itemset).
    if((k - 1) == 1):
        previous_1_itemset = prev_candidate_itemset
    # Else prev_candidate_itemset is the Lk_itemset (Large k-Itemset).
    elif((k - 1) > 1):
        
        print("need to code k >= 3")
        print(prev_candidate_itemset)
        
        # A dictionary datastructure to get the set of previous_1_itemset
        previous_1_itemset_dict = {}
        
        # Iterate through all the items in the prev_candidate_itemset.
        for prev_k_item in prev_candidate_itemset:
            prev_k_item = list(prev_k_item)
            print(prev_k_item)
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
#    print(new_candidate_k_itemset)
    del candidate_k_itemset

    return(new_candidate_k_itemset)


def prune_k_itemset(candidate_support_count,candidate_transaction_ids):
    #print(candidate_support_count)
    
    #### NEEED TO LABEL THESE DATA STRUCTURES.
    large_k_itemset = []
    large_k_support_count = {}
        
    large_k_transaction_ids = {}
    
    # Iterate over all the candidate counts.
    for k_item_key in candidate_support_count:
    
        # The support_count of the item in large_1_itemset.
        support_count = candidate_support_count[k_item_key]
        
        #print(k_item_key)
        #print(support_count)
        
        #sys.exit()
        
        # If item has a support_count less than the min_support_count then append item to item_list.
        if(support_count >= min_support_count):
            
            print(k_item_key)
            #sys.exit()
            k_item = tuple(k_item_key.split(","))
            #sys.exit()
            large_k_itemset.append(k_item)
            large_k_support_count[k_item_key] = candidate_support_count[k_item_key]
            large_k_transaction_ids[k_item_key] = candidate_transaction_ids[k_item_key]
            
    del candidate_support_count
    del candidate_transaction_ids
    
    return (large_k_itemset,large_k_support_count,large_k_transaction_ids)
    
### powerset() function

'''
Function Apriori - Apriori
'''
def apriori(transaction_database, min_support_count):

    # The frequent_itemsets dictionary data structure.
    frequent_itemsets = {}
    
    (large_1_itemset, large_1_support_count,large_1_transaction_ids) = generate_candidate_set1(transaction_database, min_support_count)

    k = 1
    
    # Generate frequent_itemsets at k = 1.
    itemset_dict = {}
    itemset_dict['itemset_data'] = large_1_itemset
    itemset_dict['itemset_support_count'] = large_1_support_count
    itemset_dict['itemset_transaction_ids'] = large_1_transaction_ids
    frequent_itemsets[str(k)] = itemset_dict
    
    print(large_1_itemset)
    print(large_1_support_count)

    
    previous_k_itemset = large_1_itemset
    k = 2
    while(len(previous_k_itemset) > 0):
        
        # Generate new apriori candidates.
        candidate_itemset_k = apriori_gen(previous_k_itemset, k)

        # The candidate in transaction counter.
        candidate_support_count = {}
        
        # The candidate in transaction ids dictionary data structure.
        candidate_transaction_ids = {}

        # Set all the k_items in candidate_itemset_k to zero in the candidate_support_count.
        for k_item in candidate_itemset_k:
            k_item = list(k_item)
            k_item_key = ",".join(k_item)
            print(k_item_key)
            candidate_support_count[k_item_key] = 0

        #print(candidate_support_count)
        #sys.exit()

        # Iterate over each transaction in the database.
        for transaction_id in transaction_database:
            print(transaction_id)
            #print(transaction_database[transaction_id])
            # Split the transaction into a list of items.
            k_itemA = transaction_database[transaction_id].split(",")
            
            # Convert the k_itemA list to a set. So we can figure out if there are subsets in candidate_itemset_k
            k_itemA = set(k_itemA)
            for k_itemB in candidate_itemset_k:
                #print(k_itemB)
                k_itemB = set(k_itemB)
                #print(k_itemB)
                # If k_itemB is a subset of k_itemA then increment
                if(k_itemB.issubset(k_itemA)):
                    print(k_itemA, k_itemB)
                    print(k_itemB.issubset(k_itemA))
                    k_itemB = sorted(list(k_itemB))
                    
                    # Join the sorted k_itemB list into a key for the candidate_support_count to increment the counter.
                    k_itemB_key = ",".join(k_itemB)
                    print(k_itemB_key)
                    
                    # Increment candidate_support_count at k_itemB_key because k_itemB_key is in the transaction.
                    candidate_support_count[k_itemB_key] = candidate_support_count[k_itemB_key] + 1
                    if(not(k_itemB_key in candidate_transaction_ids)):
                        candidate_transaction_ids[k_itemB_key] = []
                        candidate_transaction_ids[k_itemB_key].append(transaction_id)
                    elif(k_itemB_key in candidate_transaction_ids):
                        candidate_transaction_ids[k_itemB_key].append(transaction_id)
                
        print(candidate_support_count)

        print(candidate_transaction_ids)
        #sys.exit()
        
        (large_k_itemset, large_k_support_count, large_k_transaction_ids) = prune_k_itemset(candidate_support_count, candidate_transaction_ids)

        # Generate frequent_itemsets at k if large_k_itemset is not empty (length is equal to zero).
        if(len(large_k_itemset) > 0):
            itemset_dict = {}
            itemset_dict['itemset_data'] = large_k_itemset
            itemset_dict['itemset_support_count'] = large_k_support_count
            itemset_dict['itemset_transaction_ids'] = large_k_transaction_ids
            frequent_itemsets[str(k)] = itemset_dict
    
        k = k + 1
        previous_k_itemset = large_k_itemset
        
        print(frequent_itemsets)

# Minimum support count is 2.
min_support_count = 2

# Minimum confidence is 60%.
min_confidence = 0.60

# Using the "transaction_database.tsv" to test the algorithm.
transaction_database = load_transaction_database("/Users/kevin.muirhead/Desktop/macbook_air/MDSC_679/ML_Project_2/test_datasets/transaction_database1.tsv")

print(transaction_database)

apriori(transaction_database, min_support_count)

