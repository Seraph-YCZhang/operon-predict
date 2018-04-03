import numpy as np
import os
import re
from profile import work_path

def read_fna(fna_file):
    f = open(fna_file, 'r')
    fna_str = ""
    while True:
        line = f.readline().strip()
        if not line:
            break
        if ">" in line:
            pass
        else:
            fna_str = fna_str + line
    f.close()
    return fna_str


def to_fasta(gene_name, gene_str):
    gene_name = ">" + gene_name
    line_len = 80
    new_gene_str = gene_name
    i = 0
    while i < len(gene_str):
        new_gene_str = new_gene_str + "\n" + gene_str[i:i + line_len]
        i = i + line_len
    return new_gene_str


def generate_fna_according_to_gene_pos(fna_file, fna_path, gene_pos_dict):
    fna_str = read_fna(fna_file)
    fp = open(fna_path, 'w')
    fna_frg_list = []
    for key in gene_pos_dict:
        fna_frg_list.append('')
        for (start, stop) in gene_pos_dict[key]:
            fna_frg_list[-1] = fna_frg_list[-1] + fna_str[start - 1:stop]
        fp.write(to_fasta(key, fna_frg_list[-1]) + '\n')
    fp.close()


# def load_from_fasta_cds():
#     locus_tag = []
#     gene_name = []
#     f = open('eco_cds.fna','r')
#     for line in f.readlines():
#         if re.search('>',line):
#             locus_tag_tmp = re.search(r'\[locus_tag=(.{5})\]',line).group(1)
#             gene_name_tmp = re.search(r'>lcl\|(\S*)(\s)',line).group(1)
#             locus_tag.append(locus_tag_tmp)
#             gene_name.append(gene_name_tmp)
#     locus_tag = np.array(locus_tag).T
#     gene_name = np.array(gene_name).T
#     return locus_tag, gene_name

#def sort_gene(genename, gene_name, locus_tag):


def generate_fastq_list(work_path):
    fastq_list = []
    count = 0
    fastq_list_pair = ''
    for file in os.listdir(work_path):
        if re.search(r'fastq',file):
            fastq_list_pair = fastq_list_pair + str(file) + ' '
            count += 1
            if count == 4:
               fastq_list.append(fastq_list_pair.strip(' '))
               count = 0
               fastq_list_pair = ''
    return fastq_list


def kallisto_command_index(ref_seq):
    kallisto_index = ref_seq.replace('.fna','.idx')
    kallisto_index = work_path +'/' + kallisto_index
    if not os.path.exists(kallisto_index):
        os.system('kallisto index -i ' + kallisto_index + ' ' + ref_seq)
    return kallisto_index


def kallisto_command_quant(list_, _index):
    o_list=[]
    num = 0
    for pair in list_:
        tmp_list= ''
        paris = pair.split(' ')
        for file in paris:
            file = work_path + '/' + str(file)
            tmp_list = tmp_list + str(file) +' '
        o_path = work_path + '/' + 'output_' + str(num) # re.search(r'S\S{9}',pair).group(0)
        # print(o_path)
        if not os.path.exists(work_path+'/'+o_path):
            # print(tmp_list)
            command = 'kallisto quant -i ' + _index + ' -o ' +   o_path + ' ' + tmp_list + ' --pseudobam' +' > ' + str(num) + '_eco.sam'
            # print(command)
            os.system(command)
        o_list.append(o_path)
        num += 1
    return o_list


def extract_tpm_from_kallisto(work_path):
    tmp_mat = []
    for path in os.listdir(work_path):
        # print(path)
        if re.match(r'out',path):
            tmp_mat.append([])
            path_ = work_path + '/' + str(path)
            # print(path_)
            for file in os.listdir(path_):
                if re.match(r'abundance.tsv',file):
                    # print(file)
                    f = open(path_+'/'+file)
                    f_content = f.readlines()
                    for line in f_content[1:]:
                        # print(line)
                        line = line.strip()
                        line_ = line.split('\t')
                        tmp_mat[-1].append(line_[-1].strip())
                    # print(tmp_mat)
                    f.close()
    tmp_mat = np.array(tmp_mat).astype('float64').T
    return tmp_mat


def check_kallisto():
    if os.system('kallisto'):
        print('kallisto running correctly')
        return True
    else:
        print('ERROR: kallisto running incorrectly')
        return False


def check_samtools():
    if os.system('samtools'):
        print('samtools running correctly')
        return True
    else:
        print('ERROR: samtools running incorrectly')
        return False


def samtool_command(work_path):
    check = check_samtools()
    if check:
        for path in os.listdir(work_path):
            # print(path)
            m = re.match(r'(.{5}).sam',path)
            if m:
                os.system('samtools view -b '+ str(path) + ' -o ' + m.group(1) +'.bam')
                os.system('samtools sort -o' + m.group(1) + '_sort.bam ' + m.group(1) + '.bam')
                os.system('samtools index ' + m.group(1) + '_sort.bam ')
                os.system('samtools depth -a ' + m.group(1) + '_sort.bam > ' + m.group(1) + '_count.txt')


if __name__ == "__main__":
    fna_file = 'eco.fna'
    # print(load_from_fasta_cds())
    work_path = '/home/zyc/PycharmProjects/untitled'
    fastq_list = generate_fastq_list(work_path)
    ref_seq = input('the file for reference:\n')
    for file in os.listdir(work_path):
        if re.search(ref_seq, file):
            kallisto_index = kallisto_command_index(ref_seq)
    output_path = kallisto_command_quant(fastq_list, kallisto_index)
    print(extract_tpm_from_kallisto(work_path))

