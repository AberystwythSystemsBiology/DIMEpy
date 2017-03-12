import matplotlib.pyplot as plt

def roc(result, show=False, fp=False):

    y_true = result.y_true
    y_score = result.y_scores

    from sklearn.metrics import roc_curve, auc
    false_positive_ratio, true_positive_ratio, _ = roc_curve(y_true, y_score)
    area_under_curve = auc(false_positive_ratio, true_positive_ratio)
    plt.figure()
    plt.plot(false_positive_ratio, true_positive_ratio, color="k", lw=2,
             label="ROC curve (area = % 0.2f)" % area_under_curve)
    plt.xlim([0.0, 1.0])
    plt.ylim([0., 1.05])
    plt.xlabel("Specificity")
    plt.ylabel("Sensitivity")
    plt.fill_between(false_positive_ratio, true_positive_ratio, color="k", alpha=0.01)
    plt.legend(loc="lower right")
    if fp:
        plt.savefig(fp)
    if show:
        plt.show()