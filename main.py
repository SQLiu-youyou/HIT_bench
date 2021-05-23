#! /usr/bin/python2.7
# coding=utf-8

import vcf
import codecs
import sys
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np



def init_stalen():
	stalen = {'total': [0, 0, 0, 0,0,0],
            '0-99': [0, 0, 0, 0,0,0],
			'100-499': [0, 0, 0, 0,0,0],
			'500-999': [0, 0, 0, 0,0,0],
			'1k-5k': [0, 0, 0, 0,0,0],
			'5k-10k': [0, 0, 0, 0,0,0],
			'other': [0, 0, 0, 0,0,0]}
	# TP-call，call-number，TP-ans，ground-truth-number，TP-call-GT，TP-ans-GT
	return stalen

typetrans = {'insertion':'INS','deletion':'DEL','inversion':'INV','tandem':'DUP','reciprocal translocation':'BND'}
ins_dup_trans = {'INS':'DUP','DUP':'INS'}


#获取参考bed文件的gt信息
def gt_read(dictgt,path):
    f = codecs.open(path, mode='r', encoding='utf-8')
    line = f.readline()
    while line:
        a = line.split()
        line = f.readline()
        chr = a[0]
        if chr not in dictgt:
            dictgt[chr] = ''
        if float(a[-1]) > 80.0:
            dictgt[chr] = ['1/1']
        elif 80.0 >= float(a[-1]) > 20.0:
            dictgt[chr] = ['0/1']
        else:
            dictgt[chr] = ['0/0', './.']

def vcf_gt_read(dictgt,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        chr = record.CHROM
        if chr not in dictgt:
            dictgt[chr] = ''
        else:
            dictgt[chr] = record.samples[0]['GT']

#读入参考的bed文件，并建立字典dictref
def vcf_ref_read(dictref,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        sv_type = record.INFO['SVTYPE']
        chr = record.CHROM
        if sv_type not in dictref:
            dictref[sv_type] = dict()
            dictref[sv_type][chr]  = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        if sv_type == 'INS':
            dictref[sv_type][chr].append([record.POS,record.INFO['SVLEN'],0,0])
        elif sv_type == 'DEL' or sv_type == 'INV' or sv_type == 'DUP':
            dictref[sv_type][chr].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1, 0, 0])

#读入参考的bed文件，并建立字典dictref
def ref_read(dictref,path):
    f = codecs.open(path, mode='r', encoding='utf-8')
    line = f.readline()
    while line:
        a = line.split()
        # print(a[3])
        sv_type = typetrans[a[3]]
        chr = a[0]
        if sv_type not in dictref :
            dictref[sv_type] = dict()
            dictref[sv_type][chr] = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        # 字典dictref的key为染色体（str），对INS来说列表的值为：sv起始位点，变异长度,gt检测的次数,检验的次数
        if sv_type == 'INS':
            dictref[sv_type][chr].append([int(a[1]),len(a[4]),0,0])
        # 字典dictref的key为染色体（str），列表的值为：sv起始位点，终止位点以及变异长度,gt检测的次数,检验的次数
        elif sv_type=='DEL' or sv_type=='INV' or sv_type=='DUP':
            dictref[sv_type][chr].append([int(a[1]), int(a[2]),abs(int(a[2]) - int(a[1])) + 1,0,0])

        line = f.readline()

#读入检测工具call的vcf文件，并建立字典dictcall
def call_read(dictcall,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        svtype = record.INFO['SVTYPE']
        if svtype not in dictcall:
            dictcall[svtype] = dict()
            dictcall[svtype][record.CHROM] = list()
        else:
            if record.CHROM not in dictcall[svtype]:
                dictcall[svtype][record.CHROM] = list()
        #对INS来说，字典中的列表为变异起始位置，长度，GT，gt检测次数，检测次数
        if svtype == 'INS':
            dictcall[svtype][record.CHROM].append([record.POS,len(record.ALT[0]),record.samples[0]['GT'],0,0])
        #字典中的列表为变异起始位置，变异终止位置，长度，GT，gt检测次数，检测次数
        elif svtype =='DEL' or svtype== 'INV' or svtype == 'DUP':
            dictcall[svtype][record.CHROM].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1,record.samples[0]['GT'],0,0])


#对dictref和dictcall的结果进行benchmark，将所得数据填入字典的后两个域中
def benchmark(dictref,dictcall,dictgt):
    for svtype in dictcall:
        if svtype not in dictref:
            continue
        else:
            for chr in dictcall[svtype]:
                if chr not in dictref[svtype]:
                    continue
                else:
                    for call in dictcall[svtype][chr]:
                        for ref in dictref[svtype][chr]:
                            #对INS来说，经过rule来判断是否match上，成功则最后一个域加一
                            if svtype == 'INS':
                                if abs(call[0] - ref[0]) < threadshold:
                                    if float(min(call[1],ref[1]) / max(call[1], ref[1])) >= bias:
                                        call[-1]+=1
                                        ref[-1]+=1
                                        #如果gt也可以比对上，倒数第二个域也分别加一
                                        if chr in dictgt:
                                            if call[-3] in dictgt[chr]:
                                                call[-2] += 1
                                                ref [-2] += 1
                                        break
                            elif svtype =='DEL' or svtype == 'INV' or svtype == 'DUP':
                                if max(ref[0] - threadshold, call[0]) < min(ref[1] + threadshold,call[1]):
                                    if float(min(call[2], ref[2]) /
                                     max(call[2], ref[2])) >= bias:
                                        call[-1] += 1
                                        ref[-1] += 1
                                        if chr in dictgt:
                                            if call[-3] in dictgt[chr]:
                                                call[-2] += 1
                                                ref[-2] += 1
                                        break

#为下一个statistics函数服务的小函数，目的是填充svtype字典
def sv_len (lenth,n,svtype,svsize):
    if 0 < lenth < 100:
        svsize[svtype]['0-99'][n] += 1
    elif lenth < 500:
        svsize[svtype]['100-499'][n] += 1
    elif lenth < 1000:
        svsize[svtype]['500-999'][n] += 1
    elif lenth < 5000:
        svsize[svtype]['1k-5k'][n] += 1
    elif lenth < 10000:
        svsize[svtype]['5k-10k'][n] += 1
    else:
        svsize[svtype]['other'][n] += 1

#输入dictref和dictcall来统计各个类型的sv_len
def statistics(dictref,dictcall,svsize1):
    for svtype in dictcall:
        if svtype in svsize1:
            for chr in dictcall[svtype]:
                for call in dictcall[svtype][chr]:
                    #call-number的总体值
                    svsize1[svtype]['total'][1] += 1
                    #call-number各个型号的总体值
                    sv_len(call[-4],1,svtype,svsize1)
                    if call[-1] > 0:
                        #TP-call的总体值
                        svsize1[svtype]['total'][0] += 1
                        #TP-call各个型号的值
                        sv_len(call[-4],0,svtype,svsize1)
                    if call [-2] > 0:
                        #TP-call（GT）的总体值
                        svsize1[svtype]['total'][4] += 1
                        #TP-call（GT）各个型号的值
                        sv_len(call[-4], 4, svtype,svsize1)

    for svtype in dictref:
        if svtype in svsize1:
            for chr in dictref[svtype]:
                for ref in dictref[svtype][chr]:
                    # ref-number的总体值
                    svsize1[svtype]['total'][3] += 1
                    # ref-number各个型号的值
                    sv_len(ref[-3], 3, svtype,svsize1)
                    if ref[-1] > 0:
                        # TP-ref的总体值
                        svsize1[svtype]['total'][2] += 1
                        # TP-ref各个型号的值
                        sv_len(ref[-3], 2, svtype,svsize1)
                    if ref [-2] > 0:
                        # TP-ref（GT）的总体值
                        svsize1[svtype]['total'][5] += 1
                        # TP-ref（GT）各个型号的值
                        sv_len(ref[-3], 5, svtype,svsize1)


#打印precision，recall，f-score等评测值
def print_result(sv_len,F1):
    print('\t\t''pre', '\t''recall', '\t''F1', '\t''MCC')
    for key in sv_len:
        if sv_len[key][1] != 0:
            precision = float(sv_len[key][0] / sv_len[key][1])
        else:
            precision = 0
        if sv_len[key][3] != 0:
            recall = float(sv_len[key][2] / sv_len[key][3])
        else:
            recall = 0
        if precision+recall == 0:
            f_score = 0
        else:
            f_score = 2 * precision * recall / (precision + recall)
        FP = sv_len[key][1] - sv_len[key][0]
        FN = sv_len[key][3] - sv_len[key][2]
        TP = sv_len[key][0]
        TN = sv_len[key][2]
        if TP+FP == 0 or TP+FN == 0 or TN+FP == 0 or TN+FN == 0 :
            mcc = -1
        else:
            mcc = (TN*TP-FN*FP)/pow((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN),0.5)
        print('{:<15}'.format(key), '%.4f' % precision, '\t''%.4f' % recall, '\t''%.4f' % f_score,'\t''%.3f' % mcc,
              '\t''[', '{:<6}{:<6}{:<6}'.format(sv_len[key][0], sv_len[key][1], sv_len[key][2]), sv_len[key][3], ']')
        if key == 'total':
            F1.append(f_score)

#打印考虑gt的precision，recall，f-score等评测值
def print_result_gt(sv_len,F1_gt):
    print('\t\t''pre', '\t''recall', '\t''F1','\t''MCC')
    for key in sv_len:
        if sv_len[key][1] != 0:
            precision = float(sv_len[key][4] / sv_len[key][1])
        else:
            precision = 0
        if sv_len[key][3] != 0:
            recall = float(sv_len[key][5] / sv_len[key][3])
        else:
            recall = 0
        if precision+recall == 0:
            f_score = 0
        else:
            f_score = 2 * precision * recall / (precision + recall)
        FP = sv_len[key][1] - sv_len[key][4]
        FN = sv_len[key][3] - sv_len[key][5]
        TP = sv_len[key][4]
        TN = sv_len[key][5]
        if TP+FP == 0 or TP+FN == 0 or TN+FP == 0 or TN+FN == 0 :
            mcc = -1
        else:
            mcc = (TN*TP-FN*FP)/pow((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN),0.5)
        print('{:<15}'.format(key), '%.4f' % precision, '\t''%.4f' % recall, '\t''%.4f' % f_score,'\t''%.3f' % mcc,
              '\t''[', '{:<6}{:<6}{:<6}'.format(sv_len[key][4], sv_len[key][1], sv_len[key][5]), sv_len[key][3], ']')
        if key == 'total':
            F1_gt.append(f_score)

#考虑call中的dup比对到ref的ins情况
def ref_ins_call_dup():
    for svtype in dictref:
        # ref中的变异类型是INS
        if svtype == 'INS':
            for chr in dictref[svtype]:
                if chr not in dictcall['DUP']:
                    continue
                else:
                    for ref in dictref[svtype][chr]:
                        # call中的变异类型是DUP
                        for call in dictcall['DUP'][chr]:
                            #第一个筛选条件：起始位点和阈值之差，要比终止位点和阈值之和要小（对INS来说，用起始位点加变异长度来体现终止位置）
                            if max(ref[0]-threadshold,call[0])<min(ref[0]+ref[1]+threadshold,call[1]):
                                #第二个筛选条件：长度之比大于偏移量
                                if float(min(call[2], ref[1]) /max(call[2], ref[1]))>= bias:
                                    #第三个筛选条件：call中的dup未被当作dup找到过
                                    if call[-1] == 0:
                                        # print(call)
                                        call[-1]+=1
                                        ref[-1]+=1
                                        if call[-3] in dictgt[chr]:
                                            call[-2] += 1
                                            ref[-2] +=1

#考虑call中的ins比对到ref的dup情况
def ref_dup_call_ins():
    for svtype in dictref:
        # ref中的变异类型是INS
        if svtype == 'DUP':
            for chr in dictref[svtype]:
                if chr not in dictcall['INS']:
                    continue
                else:
                    for ref in dictref[svtype][chr]:
                        # call中的变异类型是INS
                        for call in dictcall['INS'][chr]:
                            #第一个筛选条件：起始位点和阈值之差，要比终止位点和阈值之和要小（对INS来说，用起始位点加变异长度来体现终止位置）
                            if max(ref[0]-threadshold,call[0]) < min(ref[1]+threadshold,call[0]+call[1]):
                                #第二个筛选条件：长度之比大于偏移量
                                if float(min(call[1],ref[2])/max(call[1],ref[2])) >= bias:
                                    #第三个筛选条件：call中的ins未被当作isn找到过
                                    if call[-1] == 0:
                                        call[-1] += 1
                                        ref[-1] += 1
                                        if call[-3] in dictgt[chr]:
                                            call[-2] += 1
                                            ref[-2] += 1

#打印ins和dup互相比对之后的各种评测值
def print_result_addition(str,F1):
    print('\t\t''pre', '\t''recall', '\t''F1','\t''MCC')
    for type in svsize1[str]:
        if svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][0] != 0 :
            precision=float(svsize1[str][type][0] + svsize2[str][type][0]) / (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][0])
        else:
            precision = 0
        if svsize1[str][type][3] != 0:
            recall=float(svsize1[str][type][2] + svsize2[str][type][2]) / (svsize1[str][type][3])
        else:
            recall = 0
        if precision+recall == 0:
            f_score = 0
        else:
            f_score = 2 * precision * recall / (precision + recall)
        FP = (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][0])-(svsize1[str][type][0] + svsize2[ins_dup_trans[str]][type][0])
        FN = svsize1[str][type][3]-(svsize1[str][type][2] + svsize2[str][type][2])
        TP = (svsize1[str][type][0] + svsize2[ins_dup_trans[str]][type][0])
        TN = (svsize1[str][type][2] + svsize2[str][type][2])
        if TP+FP == 0 or TP+FN == 0 or TN+FP == 0 or TN+FN == 0 :
            mcc = -1
        else:
            mcc = (TN*TP-FN*FP)/pow((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN),0.5)
        print('{:<15}'.format(type), '%.4f' % precision, '\t''%.4f' % recall, '\t''%.4f' % f_score,'\t''%.3f' % mcc,
              '\t''[', '{:<6}{:<6}{:<6}'.format((svsize1[str][type][0] + svsize2[ins_dup_trans[str]][type][0]), (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][0]), (svsize1[str][type][2] + svsize2[str][type][2])), svsize1[str][type][3], ']')
        if type == 'total':
            F1.append(f_score)

#打印考虑gt的ins和dup互相比对之后的各种评测值
def print_result_addition_gt(str,F1_gt):
    print('\t\t''pre', '\t''recall', '\t''F1','\t''MCC')
    for type in svsize1[str]:
        if svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][4] != 0 :
            precision=float(svsize1[str][type][4] + svsize2[str][type][4]) / (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][4])
        else:
            precision = 0
        if svsize1[str][type][3] != 0:
            recall=float(svsize1[str][type][5] + svsize2[str][type][5]) / (svsize1[str][type][3])
        else:
            recall = 0
        if precision+recall == 0:
            f_score = 0
        else:
            f_score = 2 * precision * recall / (precision + recall)
        FP = (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][4])-(svsize1[str][type][4] + svsize2[ins_dup_trans[str]][type][4])
        FN = svsize1[str][type][3] - (svsize1[str][type][5] + svsize2[str][type][3])
        TP = (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][4])
        TN = svsize1[str][type][3]
        if TP + FP == 0 or TP + FN == 0 or TN + FP == 0 or TN + FN == 0:
            mcc = -1
        else:
            mcc = (TN * TP - FN * FP) / pow((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN), 0.5)
        print('{:<15}'.format(type), '%.4f' % precision, '\t''%.4f' % recall, '\t''%.4f' % f_score,'\t''%.3f' % mcc ,
              '\t''[', '{:<6}{:<6}{:<6}'.format((svsize1[str][type][4] + svsize2[ins_dup_trans[str]][type][4]), (svsize1[str][type][1] + svsize2[ins_dup_trans[str]][type][4]), (svsize1[str][type][5] + svsize2[str][type][3])), svsize1[str][type][3], ']')
        if type == 'total':
            F1_gt.append(f_score)


# 计算基因频率 并将其写进vcf中
def Af_calculate(dictcall,path):
    try:
        vcf_reader = vcf.Reader(filename=path)
    except:
        print("file error")
        exit(0)
    AF_list = list()
    n = 0
    for record in vcf_reader:
        AF = 0
        n = len(record.samples)
        for i in range(0,len(record.samples)):
            # print(record.samples[i]['GT'])
            if record.samples[i]['GT'] == '1/1':
                AF += 2
            elif record.samples[i]['GT'] == '0/1':
                AF += 1
        AF = AF/(2 * len(record.samples))
        AF = round(AF,4)
        AF_list.append(AF)

    #写进文件
    f = codecs.open(path, mode='r', encoding='utf-8')
    line = f.readline()
    n = len(record.samples)
    i = 0
    file_handle = open(sys.argv[3], mode='w')
    INFO = 0
    while line:
        a = line.split()
        if a[0][0] == "#" and a[0][2] != "I" and a[0][2] != "F":
            file_handle.write(line)
        elif a[0][0] == "#" and a[0][2] == "I":
            file_handle.write(line)
        elif a[0][0] == "#" and a[0][2] == "F" and INFO == 0 :
            file_handle.write('##INFO=<ID=AF,Number=1,Type=Float,Description='+'"Allele frequency"'+'>'+'\n')
            file_handle.write(line)
            INFO = 1
        elif a[0][0] == "#" and a[0][2] == "F" and INFO != 0 :
            file_handle.write(line)
        elif a[0][0] != "#":
            b = a[7] + ';AF='+ str(AF_list[i])
            i += 1
            c = ""
            for x in range(0,7):
                c += a[x]
                c += "\t"
            c += b
            c += "\t"
            for y in range(8,n+8+1):
                c += a[y]
                c += "\t"
            file_handle.write(c+'\n')
        line = f.readline()

# 读入gtf文件，并为其创建字典
def	gtf_read(dictgtf,path):
    file = open(path,"r")
    for line in file:
        seq = line.strip('\n').split('\t')
        chr = seq[0][3:]
        start = seq[3]
        end = seq[4]
        infor = seq[8]
        if chr not in dictgtf:
            dictgtf[chr] = list()
        dictgtf[chr].append([start,end,infor])

# 对于INS 直接判断插入位点是否在gtf文件中即可（二分法）
def add_anotation_in(var_list, breakpoint):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = set()
    while left < right:
        mid = (left + right + 1) >> 1
        if int(var_list[mid][0]) <= breakpoint:
            left = mid
        else:
            right = mid - 1
    for i in range(left, -1, -1):
        if breakpoint < int(var_list[i][1]):
            if var_list[i][2] not in info:
                info.add(var_list[i][2])
        if breakpoint - int(var_list[i][0]) > 3000000:
            break
    return info

# 对于其他类型的SV 需要判断overlap
def add_anotation_overlap(var_list, start, end):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = set()
    while left < right:
        mid = (left + right + 1) >> 1
        if int(var_list[mid][0]) <= start:
            left = mid
        else:
            right = mid - 1
    for i in range(left, -1, -1):
        if start < int(var_list[i][1]):
            if var_list[i][2] not in info:
                info.add(var_list[i][2])
        if start - int(var_list[i][0]) > 3000000:
            break
    for i in range(left + 1, len(var_list), 1):
        if int(var_list[i][0]) < end:
            if var_list[i][2] not in info:
                info.add(var_list[i][2])
        else:
            break
    return info


def gtf_call_read(dictcallgtf,dictgtf,path):
    vcf_reader = vcf.Reader(filename=path)
    annotation = list()
    for record in vcf_reader:
        if record.INFO['SVTYPE'] =='INS':
            if record.CHROM in dictgtf:
                annotation.append(add_anotation_in(dictgtf[record.CHROM], record.POS))
            else:
                annotation.append(set())

        else:
            if record.CHROM in dictgtf:
                annotation.append(add_anotation_overlap(dictgtf[record.CHROM],record.POS,record.INFO['END']))
            else:
                annotation.append(set())
        n = len(record.samples)

    file_handle = open(sys.argv[4], mode='w')
    i = 0
    f = codecs.open(path, mode='r', encoding='utf-8')
    line = f.readline()
    while line:
        a = line.split()
        if a[0][0] == "#":
            file_handle.write(line)

        elif a[0][0] != "#":
            b = a[7] + ';'
            for x in annotation[i]:
                b += x
            i += 1
            c = ""
            for x in range(0, 7):
                c += a[x]
                c += "\t"
            c += b
            c += "\t"
            for y in range(8, n + 8 + 1):
                c += a[y]
                c += "\t"
            file_handle.write(c + '\n')
        line = f.readline()

def photo_plot(dimension,F1,F1_gt):
    bar_width = 0.3
    index_presence = np.arange(len(dimension))
    index_gt = index_presence + bar_width
    plt.bar(index_presence, height=F1, width=bar_width, color='steelblue', label='presence')
    plt.bar(index_gt, height=F1_gt, width=bar_width, color='palevioletred', label='genotype')
    # plt.legend()
    # index = np.arange(len(dimension));
    # for a, b in zip(index, F1):  # 柱子上的数字显示
    #     plt.text(a, b, '%.2f' % b, ha='center', va='bottom', fontsize=7);
    # for a, b in zip(index + bar_width, F1_gt):
    #     plt.text(a, b, '%.2f' % b, ha='center', va='bottom', fontsize=7);
    #
    plt.xticks(index_presence + bar_width / 2, dimension)
    # plt.ylabel('F1-score')
    # plt.title('F1 score under different SV types')
    plt.savefig('F1.png')

if __name__ == '__main__':

    if sys.argv[1] == 'addgt':
        # print("输入群体vcf文件名")
        afpath = sys.argv[2]
        if afpath == '-h':
            print("\n")
            print("---------------------------------usage---------------------------------")
            print("usage: python main.py addgt input.vcf output.vcf")
            print("-----------------------------------------------------------------------")
            print("\n")
        else:
            dictcallaf = dict()
            print("AF loading...")
            Af_calculate(dictcallaf, afpath)


    if sys.argv[1] == 'bench':
        # print("输入ref的GT路径")
        if sys.argv[2] == '-h':
            print("\n")
            print("---------------------------------usage---------------------------------")
            print("usage: python main.py bench region.bed base.bed comp.vcf [normal/convert]")
            print("\n")
            print("region.bed is a bed file which contains the gt information of the base.bed.")
            print("base.bed is the base calling.")
            print("comp.vcf is the comp calling.")
            print("normal indicates no conversion between insertions and duplications.")
            print("-----------------------------------------------------------------------")
            print("\n")
        else:
            gtpath = sys.argv[2]
            dictgt = dict()
            gt_read(dictgt, gtpath)

            # print("输入ans的路径")
            refpath = sys.argv[3]
            dictref = dict()
            ref_read(dictref, refpath)

            threadshold = 100  # 阈值
            bias = 0.7  # 长度偏移量

            # print("输入call的路径")
            callpath = sys.argv[4]
            dictcall = dict()
            call_read(dictcall, callpath)

            svsize1 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen()}
            benchmark(dictref, dictcall, dictgt)
            statistics(dictref, dictcall, svsize1)
            total_perform = init_stalen()
            for key in total_perform:
                for i in range(0, 6):
                    total_perform[key][i] = svsize1['INS'][key][i] + svsize1['DEL'][key][i] + \
                                            svsize1['DUP'][key][i] + svsize1['INV'][key][i]

            # print(total_perform)
            if sys.argv[5] == 'normal':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV')
                F1 = list()
                F1_gt = list()
                print("No conversion between insertions and duplications")
                print("---------------------------------TOTAL---------------------------------")
                print_result(total_perform,F1)
                print("---------------------------------INS---------------------------------")
                print_result(svsize1['INS'],F1)
                print("---------------------------------DEL---------------------------------")
                print_result(svsize1['DEL'],F1)
                print("---------------------------------DUP---------------------------------")
                print_result(svsize1['DUP'],F1)
                print("---------------------------------INV---------------------------------")
                print_result(svsize1['INV'],F1)

                print("---------------------------------TOTAL---------------------------------")
                print_result_gt(total_perform,F1_gt)
                print("--------------------------------INS-GT-------------------------------")
                print_result_gt(svsize1['INS'],F1_gt)
                print("--------------------------------DEL-GT-------------------------------")
                print_result_gt(svsize1['DEL'],F1_gt)
                print("--------------------------------DUP-GT-------------------------------")
                print_result_gt(svsize1['DUP'],F1_gt)
                print("--------------------------------INV-GT-------------------------------")
                print_result_gt(svsize1['INV'],F1_gt)

                photo_plot(dimension,F1,F1_gt)


            if sys.argv[5] == 'convert':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV')
                F1 = list()
                F1_gt = list()
                print("Converting between insertions and duplications")
                svsize2 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen()}
                ref_ins_call_dup()
                ref_dup_call_ins()

                statistics(dictref, dictcall, svsize2)
                total_perform = init_stalen()
                for key in total_perform:
                    for i in range(0, 6):
                        total_perform[key][i] = svsize2['INS'][key][i] + svsize2['DEL'][key][i] + \
                                                svsize2['DUP'][key][i] + svsize2['INV'][key][i]
                for key in svsize2:
                    for type in svsize2[key]:
                        for i in range(0, 6):
                            svsize2[key][type][i] = svsize2[key][type][i] - svsize1[key][type][i]
                print("---------------------------------TOTAL---------------------------------")
                print_result(total_perform, F1)
                print("---------------------------------INS---------------------------------")
                print_result_addition('INS',F1)
                print("---------------------------------DEL---------------------------------")
                print_result(svsize1['DEL'],F1)
                print("---------------------------------DUP---------------------------------")
                print_result_addition('DUP',F1)
                print("---------------------------------INV---------------------------------")
                print_result(svsize1['INV'],F1)

                print("---------------------------------TOTAL---------------------------------")
                print_result_gt(total_perform, F1_gt)
                print("--------------------------------INS-GT-------------------------------")
                print_result_addition_gt('INS',F1_gt)
                print("--------------------------------DEL-GT-------------------------------")
                print_result_gt(svsize1['DEL'],F1_gt)
                print("--------------------------------DUP-GT-------------------------------")
                print_result_addition_gt('DUP',F1_gt)
                print("--------------------------------INV-GT-------------------------------")
                print_result_gt(svsize1['INV'],F1_gt)

                photo_plot(dimension, F1, F1_gt)

    if sys.argv[1] == 'consistency':
        # print("输入ref的GT路径")
        if sys.argv[2] == '-h':
            print("\n")
            print("---------------------------------usage---------------------------------")
            print("usage: python main.py consistency base.vcf comp.vcf [normal/convert]")
            print("\n")
            print("base.vcf is the base calling.")
            print("comp.vcf is the comp calling.")
            print("normal indicates no conversion between insertions and duplications.")
            print("-----------------------------------------------------------------------")
            print("\n")
        else:
            gtpath = sys.argv[2]
            dictgt = dict()
            vcf_gt_read(dictgt, gtpath)

            # print("输入ans的路径")
            refpath = sys.argv[2]
            dictref = dict()
            vcf_ref_read(dictref, refpath)

            threadshold = 100  # 阈值
            bias = 0.7  # 长度偏移量

            # print("输入call的路径")
            callpath = sys.argv[3]
            dictcall = dict()
            call_read(dictcall, callpath)

            svsize1 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen()}
            benchmark(dictref, dictcall, dictgt)
            statistics(dictref, dictcall, svsize1)
            total_perform = init_stalen()
            for key in total_perform:
                for i in range(0, 6):
                    total_perform[key][i] = svsize1['INS'][key][i] + svsize1['DEL'][key][i] + \
                                            svsize1['DUP'][key][i] + svsize1['INV'][key][i]

            # print(total_perform)
            if sys.argv[4] == 'normal':
                F1 = list()
                F1_gt = list()
                print("No conversion between insertions and duplications")
                print("---------------------------------TOTAL---------------------------------")
                print_result(total_perform,F1)
                print("---------------------------------INS---------------------------------")
                print_result(svsize1['INS'],F1)
                print("---------------------------------DEL---------------------------------")
                print_result(svsize1['DEL'],F1)
                print("---------------------------------DUP---------------------------------")
                print_result(svsize1['DUP'],F1)
                print("---------------------------------INV---------------------------------")
                print_result(svsize1['INV'],F1)
                print("---------------------------------TOTAL---------------------------------")
                print_result_gt(total_perform,F1_gt)
                print("--------------------------------INS-GT-------------------------------")
                print_result_gt(svsize1['INS'],F1_gt)
                print("--------------------------------DEL-GT-------------------------------")
                print_result_gt(svsize1['DEL'],F1_gt)
                print("--------------------------------DUP-GT-------------------------------")
                print_result_gt(svsize1['DUP'],F1_gt)
                print("--------------------------------INV-GT-------------------------------")
                print_result_gt(svsize1['INV'],F1_gt)


            if sys.argv[4] == 'convert':
                F1 = list()
                F1_gt = list()
                print("Converting between insertions and duplications")
                svsize2 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen()}
                ref_ins_call_dup()
                ref_dup_call_ins()
                statistics(dictref, dictcall, svsize2)
                total_perform = init_stalen()
                for key in total_perform:
                    for i in range(0, 6):
                        total_perform[key][i] = svsize2['INS'][key][i] + svsize2['DEL'][key][i] + \
                                                svsize2['DUP'][key][i] + svsize2['INV'][key][i]
                for key in svsize2:
                    for type in svsize2[key]:
                        for i in range(0, 6):
                            svsize2[key][type][i] = svsize2[key][type][i] - svsize1[key][type][i]
                print("---------------------------------TOTAL---------------------------------")
                print_result(total_perform, F1)
                print("---------------------------------INS---------------------------------")
                print_result_addition('INS',F1)
                print("---------------------------------DEL---------------------------------")
                print_result(svsize1['DEL'],F1)
                print("---------------------------------DUP---------------------------------")
                print_result_addition('DUP',F1)
                print("---------------------------------INV---------------------------------")
                print_result(svsize1['INV'],F1)

                print("---------------------------------TOTAL---------------------------------")
                print_result_gt(total_perform, F1_gt)
                print("--------------------------------INS-GT-------------------------------")
                print_result_addition_gt('INS',F1_gt)
                print("--------------------------------DEL-GT-------------------------------")
                print_result_gt(svsize1['DEL'],F1_gt)
                print("--------------------------------DUP-GT-------------------------------")
                print_result_addition_gt('DUP',F1_gt)
                print("--------------------------------INV-GT-------------------------------")
                print_result_gt(svsize1['INV'],F1_gt)

    if sys.argv[1] == 'annotation':
        # print("输入gtf的路径")
        if sys.argv[2] == '-h':
            print("\n")
            print("---------------------------------usage---------------------------------")
            print("usage: python main.py annotation gtf input.vcf output.vcf")
            print("-----------------------------------------------------------------------")
            print("\n")
        else:
            gtfpath = sys.argv[2]
            dictgtf = dict()
            gtf_read(dictgtf, gtfpath)
            print("annotating")
            gtfcallpath = sys.argv[3]
            dictcallgtf = dict()
            gtf_call_read(dictcallgtf, dictgtf, gtfcallpath)


    if sys.argv[1] == '-h':
        print("\n")
        print("modules:bench, consistency, addgt, annotation")
        print("-------------------------------------------------Modules-------------------------------------------------")
        print("\tbench\t\t\tThe benchmark reslut will show in the monitor.")
        print("\tconsistency\t\tThis module could calculate the consistency of two VCF files. ")
        print("\taddgt\t\t\tAF can be calculate by giving a VCF.")
        print("\tannotation\t\tFor a given VCF, the annotation of every SV can be added in this module.")
        print("---------------------------------------------------------------------------------------------------------")
        print("The specific command for each moduel are available when adding '-h' to each module")
        print("This program was developed by Shiqi Liu")
        print("\n")


