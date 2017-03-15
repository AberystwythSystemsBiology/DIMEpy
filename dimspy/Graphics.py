import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec


def box_plots(da, class_df, fp="test.pdf"):
    from matplotlib.backends.backend_pdf import PdfPages
    df = da._append_class(class_df)
    classes = list(set(df[df.columns[0]].values))
    variables = df.columns[1:]
    with PdfPages(fp) as pdf_pages:
        for v_index in range(0, len(variables), 8):
            plt.figure(figsize=(9,6))
            gs = gridspec.GridSpec(3,3)
            for gs_index in range(9):
                try:
                    data = []
                    for c in classes:
                        data.append(df.ix[df[df[df.columns[0]] == c].index][variables[v_index+gs_index]].values)
                    plt.subplot(gs[gs_index])
                    plt.title(variables[v_index+gs_index], fontsize=10)
                    plt.boxplot(data, widths=0.65, labels=classes)
                    plt.xlabel(df.columns[0], fontsize=8)
                    plt.ylabel("processed intensity \n (m/z)", fontsize=8)
                    plt.yticks(fontsize=8)
                    plt.xticks(fontsize=10)
                except IndexError:
                    pass
            plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            pdf_pages.savefig()

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