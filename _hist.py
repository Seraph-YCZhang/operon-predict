import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from profile import work_path


def get_count_with_gene(position_dict, base_count):
    # gene_list = []

    count = []
    for key in position_dict:
        array_ = np.array([])
        for (x,y) in position_dict[key]:
            # print(base_count[x - 1: y])
            array_ = np.append(array_, base_count[x - 1: y])
        # gene_list.append(key)
        count.append(array_)
    return count


def count_per_gene(gene, count):
        plt.figure(dpi=100)
        plt.style.use('ggplot')
        plt.title("kernel density estimation(Gaussian):" + gene)
        plt.xlabel("count")
        plt.ylabel("f(count)")
        # print(count)
        count_ = np.unique(count)
        if count_.shape[0] > 1:
            sns.kdeplot(count, data2=None, kernel='gau', bw="scott",
                            label='scott', linewidth=0.3, shade=True)
            sns.kdeplot(count, data2=None, kernel='gau', bw="silverman",
                            label='silverman', linewidth=0.3, shade=True)
        plt.savefig(work_path + '/hist_pics/'+ str(gene))
        plt.close()
