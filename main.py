import os
import sys
import re
import numpy as np
import operator
import time
import math
import scipy.stats as stats
import matplotlib.pyplot as plt
import load_Kallisto as load_k
import _hist
from multiprocessing import Pool
from profile import fna_file, fna_new, work_path

np.set_printoptions(suppress=True)
np.seterr(divide='ignore', invalid='ignore')

def check_work_path():
    global work_path
    if os.path.exists(work_path):
        print('workpath sets rightly\n')
    else:
        print('wrong workpath')
        sys.exit(0)
    if work_path.endswith('/'):
        work_path = work_path[:-1]



def srr_download(download_choose):
    while download_choose is 'y':
        file_name = input('filename: ')
        os.system('fastq-dump --split-files --gzip '+ file_name)
        download_choose = input('next data to download?[y/n]')
        if download_choose is 'n':
            break
    return 0


def gff_download(choose):
    download_path = input('path you want to download files: ')
    print("coding...")
    return download_path


def text_read(filename):
    file = open(filename,'r')
    content = file.readlines()
    new = []
    for line in content:
        line = line.strip('\n')
        new.append(line)
    file.close()
    length = len(new)
    return new, length


def load_file(work_path):
    contentofsam, lofsam, direction, contentofgff, lofgff, base_count, genename, sposition, eposition, explvofgene, condition = [], [], [], [], [], [], [], [], [], [], []
    countfilecount, gffcount = 0, 0
    for file in os.listdir(work_path):
        if re.search(r'count',file):
            tmp1,tmp2 = text_read(file)
            contentofsam.append(tmp1)
            lofsam.append(tmp2)
            tmp3 = exract_expression_level(contentofsam[countfilecount],lofsam[countfilecount])
            base_count.append(tmp3)
            countfilecount = countfilecount+1
        elif re.search(r'gff',file):
            tmp1,tmp2 = text_read(file)
            contentofgff.append(tmp1)
            lofgff.append(tmp2)
            tmp3, tmp4, tmp5, tmp6 = extract_position(contentofgff[gffcount],lofgff[gffcount])
            genename.append(tmp3)
            sposition.append(tmp4)
            eposition.append(tmp5)
            direction.append(tmp6)
            gffcount = gffcount + 1
    return base_count,countfilecount,genename,sposition,eposition,direction,gffcount


def exract_expression_level(content, length): # from .count
    expcount = []
    for i in range(length):
        line = content[i].split()
        expcount.append(int(line[2]))
    return expcount


def extract_position(content, length): # from .gff
    startposition = []
    endposition = []
    genename = []
    midposition = []
    genelist = []
    direction = []
    directionlist = []
    count_ = 0
    for i in range(length):
        line = content[i].split()
        tmp_key = re.match(r'b(.*)', line[0])
        if tmp_key:
            if tmp_key.group(0) not in genename:
                startposition.append(line[4])
                endposition.append(line[5])
                genename.append(line[0])
                direction.append(line[8])
                midposition.append(int((startposition[count_] + endposition[count_]))/2)
                count_ = count_+1
            elif tmp_key.group(0) in genename:
                no_ = genename.index(tmp_key.group(0))
                midposition[no_] = (int(startposition[no_])+int(line[5]))/2
    print(count_)
    tmplist = zip(genename, midposition, direction)
    tmp = list(tmplist)
    sorted(tmp, key=lambda a:a[1])
    zip(tmp)
    for i in tmp:
        genelist.append(i[0])
        directionlist.append(i[2])
    return genelist, startposition, endposition, directionlist


def load_from_fasta_cds():
    locus_tag = []
    gene_name = []
    f = open('eco_cds.fna','r')
    for line in f.readlines():
        if re.search('>'):
            locus_tag_tmp = re.search(r'[locus_tag=(b.*)]',line).group(1)
            gene_name_tmp = re.search(r'>lcl\|(.*)(\s){1}').group(1)
            locus_tag.append(locus_tag_tmp)
            gene_name.append(gene_name_tmp)
    return locus_tag, gene_name


def samtools_command(work_path):
    load_k.check_samtools()
    load_k.samtool_command(work_path)


def kallisto_command():
    if load_k.check_kallisto() is True:
        fastq_list = load_k.generate_fastq_list(work_path)
        load_k.generate_new_fna(fna_file, fna_new, positiondict)
        for file in os.listdir(work_path):
            if re.search(fna_new, file):
                kallisto_index = load_k.kallisto_command_index(fna_new)
        output_path = load_k.kallisto_command_quant(fastq_list, kallisto_index)
        mat = load_k.extract_tpm_from_kallisto(work_path)
    return mat


def sort_gene(positiondict):
    mid_pos = []
    dict = []
    for key in positiondict:
        dict.append(key)
        num = 0
        for (x, y) in positiondict[key]:
            x_ = int(x)
            y_ = int(y)
            if num == 0:
                start = x_
                end = y_
            else:
                end = y_
            num += 1
        mid_pos.append((start + end)/2)
    dict = [dict for (mid_pos, dict) in sorted(zip(mid_pos, dict))]
    return dict


def dict_of_gene(genename, startposition, endposition):
    dict = {}
    i = 0
    for name_ in genename:
        if name_ in dict:
            dict[name_].append((int(startposition[i]), int(endposition[i])))
        else:
            dict[name_] = [(int(startposition[i]),int(endposition[i]))]
        i = i + 1
    return dict


def exp_of_gene_tpm(dict, base_count):
    explv = {}
    sum_2 = 0
    for key in dict:
        sum_1 = 0
        length = 0
        for (x,y) in dict[key]:
            sum_1 = sum_1 + sum(base_count[x-1: y])
            length = (int(length) + int(y) - int(x)) / 1000
        explv[key] = sum_1/length
        sum_2 = sum_2 + explv[key]
    for key in dict:
        explv[key] = (explv[key]/sum_2)*1000000
    explv = sorted(explv.items(),key=operator.itemgetter(0))
    return explv


def mat_choose(choose, matrix):
    if choose is '1':
        result_mat = matrix_cij(matrix) * 2
    elif choose is '2':
        result_mat = matrix_pearson(matrix)
    elif choose is '3':
        result_mat = matrix_spearman(matrix)
    else:
        print('\n!WRONG INPUT CHARACTER!')
        sys.exit(0)
    return result_mat


def pearson(vector1, vector2):
    n = vector1.shape[1]
    #simple sums
    sum1 = vector1.sum()
    sum2 = vector2.sum()
    # sum1 = sum(float(vector1[i]) for i in range(n))
    # sum2 = sum(float(vector2[i]) for i in range(n))
    #sum up the squares
    sum1_pow = (vector1**2).sum()
    sum2_pow = (vector2**2).sum()
    #sum up the products
    p_sum = (vector1*vector2).sum()
    num = p_sum - (sum1*sum2/n)
    den = math.sqrt((sum1_pow-pow(sum1, 2)/n)*(sum2_pow-pow(sum2, 2)/n))
    if den == 0:
        return 0.0
    return num/den


def matrix_spearman(matrix):
    '''
    row_ = matrix.shape[0]
    spearman_mat = []
    for i in range(row_):
        spearman_mat.append([])
        for j in range(row_):
            vec1, vec2 = matrix[i:i + 1, :], matrix[j:j + 1, :]
            #print(vec1[0])
            time_1=time.time()
            spearman_coefficient, pval = stats.spearmanr(vec1[0], vec2[0])
            print(time.time()-time_1)
            spearman_mat[i].append(spearman_coefficient)
            #print(spearman_coefficient)
    '''
    spearman_mat, p_value= stats.spearmanr(matrix, axis=1)
    return spearman_mat


def matrix_pearson(matrix):
    # row_ = matrix.shape[0]
    # pearson_mat = []
    # for i in range(row_):
    #     pearson_mat.append([])
    #     for j in range(row_):
    #         vec1, vec2 = matrix[i:i + 1, :], matrix[j:j + 1, :]
    #         pearson_coefficient = pearson(vec1, vec2)
    #         pearson_mat[i].append(pearson_coefficient)
    pearson_mat = np.corrcoef(matrix, rowvar=True)
    return pearson_mat


def matrix_cij(matrix):
    row_ = matrix.shape[0]
    column_ = matrix.shape[1]
    for i in range(column_):
        exp_i = matrix[:, i:i + 1]
        diff_exp_i = matrix - exp_i
        if i == 0:
            diff_exp = diff_exp_i
        else:
            diff_exp = np.hstack((diff_exp,diff_exp_i))
    diff_exp[diff_exp >= 0] = 1
    diff_exp[diff_exp < 0] = -1
    i_j_coexp_mat=np.empty((row_,row_))
    for i in range(row_):
        exp_i = diff_exp[i:i+1, :]
        product_of_i = diff_exp * exp_i
        mat_tmp=product_of_i.sum(axis=1)
        i_j_coexp_mat[i:i+1,:] = mat_tmp
    i_j_coexp_mat = i_j_coexp_mat/(2*column_**2)
    return i_j_coexp_mat

    # diff_exp
    #     c2-c1 c3-c1
    # g1
    # g2



def partition_with_direction(direction):
    groups=[]
    array_dirct=[]
    array_1=[]
    start=0
    i=0
    f=open('part.txt','w')
    while i < len(direction):
        array_dirct.append(direction[i])
        i += 1
        if i - start > 1:
            if array_dirct[-1] != array_dirct[0]:
                array_1.append(array_dirct[-1])
                if len(array_1) > 1 :
                    end = i - 3
                    groups.append((start,end))
                    f.write(str(direction[start:end+1]))
                    array_dirct = []
                    array_1 = []
                    start = i - 2
                    i = start
            elif array_dirct[-1] == array_dirct[0] and array_dirct[-1]!=array_dirct[-2]:
                array_1 = []
        if i == len(direction)-1:
            groups.append((start, i-1))
            f.write(str(direction[start:i + 1]))
    f.close()
    return groups


def find_max(dist_mat, threashold):
    # tmp = threashold
    x_ = 999
    max_ = max(dist_mat)
    if max_ >= threashold:
        for i in range(len(dist_mat)):
            if dist_mat[i] >= max_:
                tmp = dist_mat[i]
                x_ = i
                return x_, tmp
    if x_ == 999:
        return False


def find_nearest(x, y, matrix):
    dist_mat = []
    dist = []
    len_ = y - x
    for i in range(len_):
        dist_mat.append(matrix[x + i][x + i + 1])
        dist.append((i, i + 1))
        # if i == len_ - 1:
        #     break
    return dist_mat, dist


def check_item(x, items):
    i = 0
    for item in items:
        # print(x)
        # print(item)
        # print(items)
        if x in item:
            i += 1
            return True
    if i == 0:
        return False


def cluster(matrix, groups, genename, threashold):
    clusters = []
    for (x, y) in groups:
        check_list = []
        len_ = y - x
        gene_list = genename[x:y+1]
        cluster = []
        # print(dist)
        dist_mat, dist = find_nearest(x, y, matrix)
        print('dist: '+str(dist))
        print('dist_mat: '+str(dist_mat))
        len_2 = 0
        while len_2 < len_:
            # print(find_max(dist_mat,threashold))
            if find_max(dist_mat,threashold) is not False:
                a1, b1 = find_max(dist_mat,threashold)
                nearest_gene = dist[a1][1]
                gene_self = dist[a1][0]
                dist_mat.remove(dist_mat[a1])
                dist.remove(dist[a1])
                [gene_1,gene_2] = [gene_list[gene_self], gene_list[nearest_gene]]
                item_ = [gene_1,gene_2]
                if len_2 == 0:
                    cluster.append(item_)
                    check_list.append(item_)
                else:
                    check_gene_1 = check_item(gene_1, check_list)
                    check_gene_2 = check_item(gene_2, check_list)
                    # print('listcheck：' + str(check_list))
                    # print(check_gene_1)
                    # print(check_gene_2)
                    if check_gene_2 is True and check_gene_1 is True:
                        # print('clustercheck222：' + str(cluster))
                        for item in cluster:
                            # print('itemcheck222：' + str(item))
                            check_gene_1_initem = check_item(gene_1, item)
                            check_gene_2_initem = check_item(gene_2, item)
                            if check_gene_1_initem is True:

                                item_1 = item
                            if check_gene_2_initem is True:
                                item_2 = item
                        item__ = []
                        item__.extend(item_1)
                        item__.extend(item_2)
                        # print(cluster)
                        # print(item_1)
                        # print(item_2)
                        # print(item__)
                        cluster.remove(item_1)
                        cluster.remove(item_2)
                        # print('clustercheck：'+str(cluster))
                        cluster.append(item__)
                        check_list.append(item__)
                    elif check_gene_1 is False and check_gene_2 is False:
                        # print('1111:'+str(gene_1))
                        # print('2222:'+str(gene_2))
                        # print(cluster)
                        cluster.append(item_)
                        check_list.append(item_)
                    else:
                        for item in cluster:
                            # print('itemcheck2222111：' + str(item))
                            # print('clustercheck222：' + str(cluster))
                            # print('listcheck222：' + str(check_list))
                            check_gene_1_initem = check_item(gene_1, item)
                            check_gene_2_initem = check_item(gene_2, item)
                            if check_gene_1 is True and check_gene_2 is False and check_gene_1_initem is True:
                                # print('clustercheck221：' + str(cluster))
                                # print('listcheck221：' + str(check_list))
                                # print(check_gene_1_initem)
                                    # print(gene_2)
                                tmp = item
                                    # print('tmp:'+str(tmp))
                                tmp.append(gene_2)
                                item__ = tmp
                                # print('tmp:' + str(tmp))
                                #     # print('itemcheck：' + str(cluster))
                                # print('__:'+str(item__))
                                # print('item:'+str(item))
                                cluster.remove(item)
                                cluster.append(item__)
                                check_list.append(item__)
                                break
                                    # print('clustercheck221111：' + str(cluster))
                                    # print('listcheck221111：' + str(check_list))
                                    # check_list.append(check_gene_2)
                            if check_gene_2 is True and check_gene_1 is False and check_gene_2_initem is True:
                                # print('clustercheck2222：' + str(cluster))
                                # print('listcheck2222：' + str(check_list))
                                tmp = item
                                tmp.append(gene_1)
                                item__ = tmp
                                    # print('clustercheck：' + str(cluster))
                                    # print(item__)
                                    # print(item)
                                cluster.remove(item)
                                cluster.append(item__)
                                check_list.append(item__)
                                break
                                    # print('clustercheck：' + str(cluster))
                                    # print('listcheck：' + str(check_list))
                                    # check_list.append(check_gene_1)
            elif find_max(dist_mat, threashold) is False:
                # break
                len__=len(dist_mat)
                i = len__
                for i in range(len__):
                    nearest_gene = dist[i][1]
                    gene_self = dist[i][0]
                    i -= 1
                    gene_1, gene_2 = gene_list[gene_self], gene_list[nearest_gene]
                    check_gene_1 = check_item(gene_1, check_list)
                    check_gene_2 = check_item(gene_2, check_list)
                    # print(check_gene_1)
                    # print(check_gene_2)
                    # print(check_list)
                    # print(cluster)
                    if check_gene_1 is False and check_gene_2 is False:
                        cluster.append(gene_1)
                        cluster.append(gene_2)
                        check_list.append(gene_1)
                        check_list.append(gene_2)
                    if check_gene_1 is True and check_gene_2 is False:
                        cluster.append(gene_2)
                        check_list.append(gene_2)
                    if check_gene_2 is True and check_gene_1 is False:
                        cluster.append(gene_1)
                        check_list.append(gene_1)
            len_2 += 1
        print('cluster: '+ str(cluster))
        for item in cluster:
            clusters.append(str(item).replace('(','').replace(')','').replace(' ','').replace('\'','').replace(',',';').replace('[','').replace(']',''))
    return np.array(clusters)


def possible_cluster(partition, genename):
    list_ = []
    for (x, y) in partition:
        len_ = y - x
        count = 0
        while count < len_ + 1:
            for j in range(count + 1):
                len__ = len_ - count
                item_ = ''
                for n in range(j, len__ + 1):
                    item_=item_+str(genename[x + n].strip('\n'))+';'
                list_.append(str(item_).replace(',',';').replace(' ','').replace('\'','').strip().strip(';'))
            count += 1
    f= open('ps.txt','w+')
    for line in list_:
        f.write(line+'\n')
    f.close()
    return list_


def negative_generate(partition, genename, p_data):
    possible_data = possible_cluster(partition, genename)
    ret = [item for item in possible_data if item not in p_data]
    len_ = len(ret)
    return ret, len_


def rate_calculation(p_len, c_len, n_len, result_of_cluster, p_data, n_data):
    p_count = 0
    for item in result_of_cluster:
        if item in p_data:
            p_count += 1
        # elif item in n_data:
        #     n_count += 1
    # print(str(p_count)+' '+str(n_count)+' '+str(data_p_len))
    tpr = p_count / p_len # recall
    precison = p_count / c_len
    fpr = (c_len - p_count) / 20000
    return tpr, precison, fpr


def  roc_function(rate_list):
    plt.figure(1)
    plt.style.use('ggplot')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('FPR( False Positive Rate)')
    plt.ylabel('TPR( True Positive Rate)')
    roc_auc = 0
    prev_x = 0
    for [x, y] in rate_list:
        if x != prev_x:
            roc_auc += (x - prev_x) * y
            prev_x = x
    list_ = np.array(rate_list)
    x_ = list_[:, 0]
    # print(x_)
    y_ = list_[:, 1]
    # print(y_)
    plt.plot(x_, y_, color='r')
    # plt.savefig(threshold+'.jpg')
    plt.title('ROC Curve (AUC = %0.4f)' % roc_auc)
    plt.savefig('roc_fig.png')
    plt.show()


def  pr_function(rate_list):
    plt.figure(1)
    plt.style.use('ggplot')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    pr_auc = 0
    prev_x = 0
    for [x, y] in rate_list:
        if x != prev_x:
            pr_auc += (x - prev_x) * y
            prev_x = x
    list_ = np.array(rate_list)
    x_ = list_[:, 0]
    # print(x_)
    y_ = list_[:, 1]
    # print(y_)
    plt.plot(x_, y_, color='m')
    # plt.savefig(threshold+'.jpg')
    plt.title('PR (AUC = %0.4f)' % pr_auc)
    plt.savefig('pr_fig.png')
    plt.show()


if __name__ == "__main__":
    check_work_path()
    time_start = time.time()
    print('---operon predictor is running...')
    download_srr_choose = input('do you need to download SRR seq from ncbi[y/n]')
    if download_srr_choose is 'y':
        srr_download(download_srr_choose)
    download_gff_choose = input('do you need to download gff file from ncbi[y/n]')
    if download_gff_choose is 'y':
        gff_download(download_gff_choose)
    func_choose = input('function choose: 1 for hist_plot; 2 for operon predicion')
    if func_choose is not '1' and func_choose is not '2':
        print('WRONG INPUT!')
        sys.exit(1)
    print('---calculating...')
    #contentofsam, lofsam, direction, contentofgff, lofgff, base_count, genename, sposition, eposition, explvofgene, condition = [], [], [], [], [], [], [], [], [], [], []
    #countfilecount, gffcount = 0, 0
    base_count,countfilecount,genename,sposition,eposition,direction,gffcount = load_file(work_path)
    if gffcount == 0 :
        print('please check the files you have download')
        sys.exit(0)
    print(len(genename[0]))
    f = open ('genename.txt','w')
    for item in genename[0]:
        f.write(item+'\n')
    f.close()
    positiondict = dict_of_gene(genename[0], sposition[0], eposition[0])
    dict_gene = sort_gene(positiondict)
    # positiondict = dict_of_gene(dict_gene,)
    print(len(positiondict))
    if func_choose is '1':
        samtools_command(work_path)
        if countfilecount == 0:
            print('ERROR: u need to run kallsito function first')
            sys.exit(0)
        choose_num = str(input('no. of condition_seq: 0 - '+ str(countfilecount)))
        if not os.path.exists(work_path + '/hist_pics'):
            os.system('mkdir ' + work_path + '/hist_pics')
        # print(base_count)
        _count = _hist.get_count_with_gene(positiondict, base_count[choose_num])
        pool = Pool(10)
        result = pool.starmap(_hist.count_per_gene, zip(genename[0], _count))
        pool.close()
        pool.join()
    # my_tpm_calculation
    # for j in range(countfilecount):
    #     explvofgene.append(exp_of_gene_tpm(positiondict, base_count[j]))
    #     #print(explvofgene)
    # for i in range(countfilecount):
    #     condition.append([])
    #     for j in range(len(positiondict)):
    #         a= '%.3f' % explvofgene[i][j][1]
    #         condition[i].append(float(a))
    # print(condition)
    # matrix = np.array(condition).T
    # print(matrix)
    if func_choose is '2':
        matrix = kallisto_command()
        cor_choose = input('1 for cij; 2 for pearson; 3 for spearman\n')
        print('---calculating...')
        result_mat = mat_choose(cor_choose, matrix)
        result_of_partition = partition_with_direction(direction[0])
        # print(direction[0])
        # print(len(genename[0]))
        # print(result_of_partition)
        # print(result_of_partition)
        # f=open('dirct.txt','w')
        # f.write(str(direction[0]))
        # f.close()
        np.savetxt('result_mat.txt', result_mat, fmt="%.4f", delimiter="\t")
        positive_data, data_p_len = text_read('positiveData.txt')
        # print(positive_data)
        negative_data, data_n_len = negative_generate(result_of_partition, genename[0], positive_data)
        # print(data_n_len)
        threashold_ = input('please input the threshold(exp.: 0 1):\n')
        threashold_min = float(threashold_.split(' ')[0])
        threashold_max = float(threashold_.split(' ')[1])
        threashold = threashold_max
        minus_per_round = input('please key in the difference between ech round:(exp.: 0.05 or 0.1)\n')
        minus_per_round = float(minus_per_round)
        rate_pr_list = []
        rate_roc_list = []
        num = 0
        while threashold >= threashold_min :
            print('threshold:'+str(threashold))
            result_of_cluster = cluster(result_mat, result_of_partition, genename[0], threashold)
            f = open('result_pre'+ str(num) + '.txt','w')
            for i in result_of_cluster:
                f.write(i+'\n')
            f.close()
            c_len = len(result_of_cluster)
            print('clength:'+str(c_len))
            tpr, precision, fpr = rate_calculation(data_p_len, c_len, data_n_len, result_of_cluster, positive_data, negative_data)
            rate_pr_list.append([tpr, precision])
            rate_roc_list.append([fpr, tpr])
            threashold -= minus_per_round
            threashold = round(threashold, 2)
            num += 1
            # print(rate_list)
        pr_function(rate_pr_list)
        roc_function(rate_roc_list)
    time_end = time.time()
    print('totally cost: ', time_end - time_start)
    print('---program is over...')
