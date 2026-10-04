import numpy as np

def pairwise_auc(y, p, w=None):
    y = np.asarray(y)
    p = np.asarray(p)
    w = np.ones(len(y)) if w is None else np.asarray(w)
    pp = p[y == 1]
    pn = p[y == 0]
    wp = w[y == 1]
    wn = w[y == 0]
    comparison = (pp[:, None] > pn[None, :]).astype(float) + 0.5 * (pp[:, None] == pn[None, :])
    return np.sum(comparison * wp[:, None] * wn[None, :]) / (wp.sum() * wn.sum())

def threshold_ap(y, p, w=None):
    y = np.asarray(y)
    p = np.asarray(p)
    w = np.ones(len(y)) if w is None else np.asarray(w)
    total_positive = w[y == 1].sum()
    last_recall = 0.0
    value = 0.0
    for threshold in sorted(set(p), reverse=True):
        selected = p >= threshold
        tp = w[selected & (y == 1)].sum()
        precision = tp / w[selected].sum()
        recall = tp / total_positive
        value += precision * (recall - last_recall)
        last_recall = recall
    return value
