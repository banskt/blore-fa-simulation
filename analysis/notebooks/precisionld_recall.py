import numpy as np

def precisionld_recall_curve(ldresult, ldcut):
    #epsilon = -100
    nitems = len(ldresult)
    estimated_values_rounded = np.zeros(nitems)
    for i in range(nitems):
    #    if ldresult[i].stat >= epsilon:
        estimated_values_rounded[i] = ldresult[i].stat
    sorted_index = np.argsort(estimated_values_rounded)[::-1]
    
    real_values = np.array([x.causality for x in ldresult])
    positives = np.sum(real_values)
    negatives = nitems - positives
    
    fpr = []
    tpr = []
    precision = []
    recall = []
    nselected = []
    false_positives = 0
    true_positives  = 0
    alpha = 0
    precision_ld = []
    true_positives_with_ld = 0
    tpld_prev = 0
    fdr = []
    j = 0
    while j < nitems:
        if not estimated_values_rounded[sorted_index[j]] == alpha:
            fpr.append(false_positives / negatives)
            tpr.append(true_positives / positives)
            if (true_positives + false_positives) == 0:
                precision.append(1.0)
                recall.append(0.0)
                precision_ld.append(1.0)               
            else:
                precision.append(true_positives / (true_positives + false_positives))
                recall.append(true_positives / positives)
                precision_ld.append(true_positives_with_ld / (true_positives + false_positives))
                
            if (false_positives == 0):
                fdr.append(0.0)
            else:
                fdr.append((true_positives_with_ld - true_positives) / false_positives)
            nselected.append(j)
            alpha = estimated_values_rounded[sorted_index[j]]

        if real_values[sorted_index[j]] == 1:
            true_positives += 1
        else:
            false_positives += 1
            
        if ldresult[sorted_index[j]].ld > ldcut:
            true_positives_with_ld += 1

        j += 1

    fpr.append(false_positives / negatives)
    tpr.append(true_positives / positives)
    precision.append(true_positives / (true_positives + false_positives))
    recall.append(true_positives / positives)
    precision_ld.append(true_positives_with_ld / (true_positives + false_positives))
    fdr.append((true_positives_with_ld - true_positives) / false_positives)
    nselected.append(nitems)

    return fpr, tpr, precision, recall, nselected, precision_ld, fdr
