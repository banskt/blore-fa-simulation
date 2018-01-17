import numpy as np

def trapezoid_area(x1, x2, y1, y2):
    base = abs(x1 - x2)
    height = (y1 + y2) / 2.0
    return base * height

def roc_values(real_values, estimated_values, epsilon):
    nitems = len(estimated_values)
    estimated_values_rounded = np.zeros(nitems)
    for i in range(nitems):
        if estimated_values[i] >= epsilon:
            estimated_values_rounded[i] = estimated_values[i]
    sorted_index = np.argsort(estimated_values_rounded)[::-1]
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
    area = 0
    fp_prev = 0
    tp_prev = 0
    j = 0
    while j < nitems:
        if not estimated_values_rounded[sorted_index[j]] == alpha:
            fpr.append(false_positives / negatives)
            tpr.append(true_positives / positives)
            if (true_positives + false_positives) == 0:
                precision.append(1.0)
                recall.append(0.0)
            else:
                precision.append(true_positives / (true_positives + false_positives))
                recall.append(true_positives / positives)
            nselected.append(j)
            area += trapezoid_area(false_positives, fp_prev, true_positives, tp_prev)
            alpha = estimated_values_rounded[sorted_index[j]]
            fp_prev = false_positives
            tp_prev = true_positives

        if real_values[sorted_index[j]] == 1:
            true_positives += 1
        else:
            false_positives += 1

        j += 1

    fpr.append(false_positives / negatives)
    tpr.append(true_positives / positives)
    precision.append(true_positives / (true_positives + false_positives))
    recall.append(true_positives / positives)
    nselected.append(nitems)
    area += trapezoid_area(negatives, fp_prev, negatives, tp_prev)
    area /= (positives * negatives)

    return fpr, tpr, precision, recall, nselected, area
