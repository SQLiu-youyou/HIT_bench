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
threadshold = {'INS':800, 'DEL':800,'INV':1000,'DUP':1200,'BND':1000}

def phase_bnd(alt, chr, pos):
	# print(alt)
	# print(str(alt[0]))
	if alt[0] == ']':
		form = ']]N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	elif alt[0] == '[':
		form = '[[N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	else:
		# print(type(alt))
		if alt[1] == ']':
			form = 'N]]'
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])
		else:
			form = 'N[['
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])

	try:
		if int(chr) <= int(chr2):
			if form == 'N[[':
				form = ']]N'
			if form == ']]N':
				form = 'N[['
			return form, chr, pos, chr2, pos2
		else:
			return form, chr2, pos2, chr, pos
	except:
		return form, chr, pos, chr2, pos2


#读入带gtf的参考的vcf文件，并建立字典dictref
def gtf_vcf_ref_read(dictref,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        sv_type = record.INFO['SVTYPE']
        chr = record.CHROM
        try:
            genotype = record.samples[0]['GT']
        except:
            continue
        try:
            annotype = record.INFO['ANNOTYPE']
            # print(annotype)
        except:
            continue
        if sv_type not in dictref:
            dictref[sv_type] = dict()
            dictref[sv_type][chr]  = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        if sv_type == 'INS':
            dictref[sv_type][chr].append([record.POS,record.INFO['SVLEN'],annotype,genotype, 0,0])
        elif sv_type == 'DEL' or sv_type == 'INV' or sv_type == 'DUP':
            dictref[sv_type][chr].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1,annotype, genotype,0 , 0])
        elif sv_type == "BND":
            form, chr, start, chr2, pos2 = phase_bnd(str(record.ALT[0]), chr, record.POS)
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
            dictref[sv_type][chr].append([start, chr2, pos2, form ,annotype,genotype, 0 , 0])


def vcf_ref_read(dictref,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        sv_type = record.INFO['SVTYPE']
        chr = record.CHROM
        try:
            genotype = record.samples[0]['GT']
        except:
            continue
        if sv_type not in dictref:
            dictref[sv_type] = dict()
            dictref[sv_type][chr]  = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        if sv_type == 'INS':
            dictref[sv_type][chr].append([record.POS,record.INFO['SVLEN'],genotype, 0,0])
        elif sv_type == 'DEL' or sv_type == 'INV' or sv_type == 'DUP':
            dictref[sv_type][chr].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1, genotype,0 , 0])
        elif sv_type == "BND":
            form, chr, start, chr2, pos2 = phase_bnd(str(record.ALT[0]), chr, record.POS)
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
            dictref[sv_type][chr].append([start, chr2, pos2, form ,genotype, 0 , 0])

#对dictref和dictcall的结果进行benchmark，将所得数据填入字典的后两个域中
def benchmark(dictref,dictcall):
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
                                if abs(call[0] - ref[0]) < threadshold[svtype]:
                                    if float(min(call[1],ref[1]) / max(call[1], ref[1])) >= bias and flag==1 and anno_type in ref[-4]:
                                        call[-1]+=1
                                        ref[-1]+=1
                                        #如果gt也可以比对上，倒数第二个域也分别加一
                                        if call[-3] == ref[-3]:
                                            call[-2] += 1
                                            ref[-2] += 1
                                        break
                                    if float(min(call[1],ref[1]) / max(call[1], ref[1])) >= bias and flag==0:
                                        call[-1]+=1
                                        ref[-1]+=1
                                        #如果gt也可以比对上，倒数第二个域也分别加一
                                        if call[-3] == ref[-3]:
                                            call[-2] += 1
                                            ref[-2] += 1
                                        break
                            elif svtype =='DEL' or svtype == 'INV' or svtype == 'DUP':
                                if max(ref[0] - threadshold[svtype], call[0]) < min(ref[1] + threadshold[svtype],call[1]):
                                    if float(min(call[2], ref[2]) /
                                     max(call[2], ref[2])) >= bias and flag==1 and anno_type in ref[-4]:
                                        call[-1] += 1
                                        ref[-1] += 1
                                        if call[-3] == ref[-3]:
                                            call[-2] += 1
                                            ref[-2] += 1
                                        break
                                    if float(min(call[2], ref[2]) /
                                     max(call[2], ref[2])) >= bias and flag==0:
                                        call[-1] += 1
                                        ref[-1] += 1
                                        if call[-3] == ref[-3]:
                                            call[-2] += 1
                                            ref[-2] += 1
                                        break
                            elif svtype == 'BND' :
                                if call[1] == ref[1] and abs(call[0] - ref[0]) < threadshold[svtype] and abs(call[2] - ref[2]) < threadshold[svtype] and call[3] == ref[3] and flag==1 and anno_type in ref[-4]:
                                    call[-1] += 1
                                    ref[-1] += 1
                                    if call[-3] == ref[-3]:
                                        call[-2] += 1
                                        ref[-2] += 1
                                    break
                                if call[1] == ref[1] and abs(call[0] - ref[0]) < threadshold[svtype] and abs(call[2] - ref[2]) < threadshold[svtype] and call[3] == ref[3] and flag == 0:
                                    call[-1] += 1
                                    ref[-1] += 1
                                    if call[-3] == ref[-3]:
                                        call[-2] += 1
                                        ref[-2] += 1
                                    break
    if flag==1:
        print("OKOKOKOKOOKOOO")

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
def statistics(dictref,dictcall,svsize1,n):
    for svtype in dictcall:
        if svtype in svsize1:
            if svtype != 'BND':
                for chr in dictcall[svtype]:
                    for call in dictcall[svtype][chr]:
                        # call-number的总体值
                        svsize1[svtype]['total'][1] += 1
                        # call-number各个型号的总体值
                        sv_len(call[-4], 1, svtype, svsize1)
                        if call[-1] > 0:
                            # TP-call的总体值
                            svsize1[svtype]['total'][0] += 1
                            # TP-call各个型号的值
                            sv_len(call[-4], 0, svtype, svsize1)
                        if call[-2] > 0:
                            # TP-call（GT）的总体值
                            svsize1[svtype]['total'][4] += 1
                            # TP-call（GT）各个型号的值
                            sv_len(call[-4], 4, svtype, svsize1)
            else:
                for chr in dictcall[svtype]:
                    for call in dictcall[svtype][chr]:
                        # call-number的总体值
                        svsize1[svtype]['total'][1] += 1
                        if call[-1] > 0:
                            # TP-call的总体值
                            svsize1[svtype]['total'][0] += 1
                        if call[-2] > 0:
                            # TP-call（GT）的总体值
                            svsize1[svtype]['total'][4] += 1

    for svtype in dictref:
        if svtype in svsize1:
            if svtype!='BND':
                for chr in dictref[svtype]:
                    for ref in dictref[svtype][chr]:
                        # ref-number的总体值
                        svsize1[svtype]['total'][3] += 1
                        # ref-number各个型号的值
                        sv_len(ref[n], 3, svtype, svsize1)
                        if ref[-1] > 0:
                            # TP-ref的总体值
                            svsize1[svtype]['total'][2] += 1
                            # TP-ref各个型号的值
                            sv_len(ref[n], 2, svtype, svsize1)
                        if ref[-2] > 0:
                            # TP-ref（GT）的总体值
                            svsize1[svtype]['total'][5] += 1
                            # TP-ref（GT）各个型号的值
                            sv_len(ref[n], 5, svtype, svsize1)
            else:
                for chr in dictref[svtype]:
                    for ref in dictref[svtype][chr]:
                        # ref-number的总体值
                        svsize1[svtype]['total'][3] += 1
                        if ref[-1] > 0:
                            # TP-ref的总体值
                            svsize1[svtype]['total'][2] += 1
                        if ref[-2] > 0:
                            # TP-ref（GT）的总体值
                            svsize1[svtype]['total'][5] += 1

#打印precision，recall，f-score, mcc等评测值
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
        indicator=list()
        indicator.append(precision)
        indicator.append(recall)
        indicator.append(f_score)
        indicator.append(mcc)
        if key == 'total':
            F1.append(indicator[i])
        size_fig.append(indicator[i])

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
        indicator = list()
        indicator.append(precision)
        indicator.append(recall)
        indicator.append(f_score)
        indicator.append(mcc)
        if key == 'total':
            F1_gt.append(indicator[i])
        size_fig_gt.append(indicator[i])

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
                            if max(ref[0]-threadshold[svtype],call[0])<min(ref[0]+ref[1]+threadshold[svtype],call[1]):
                                #第二个筛选条件：长度之比大于偏移量
                                if float(min(call[2], ref[1]) /max(call[2], ref[1]))>= bias:
                                    #第三个筛选条件：call中的dup未被当作dup找到过
                                    if call[-1] == 0:
                                        # print(call)
                                        call[-1]+=1
                                        ref[-1]+=1
                                        if call[-3]==ref[3]:
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
                            if max(ref[0]-threadshold[svtype],call[0]) < min(ref[1]+threadshold[svtype],call[0]+call[1]):
                                #第二个筛选条件：长度之比大于偏移量
                                if float(min(call[1],ref[2])/max(call[1],ref[2])) >= bias:
                                    #第三个筛选条件：call中的ins未被当作isn找到过
                                    if call[-1] == 0:
                                        call[-1] += 1
                                        ref[-1] += 1
                                        if call[-3] ==ref[3]:
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
        indicator = list()
        indicator.append(precision)
        indicator.append(recall)
        indicator.append(f_score)
        indicator.append(mcc)
        if type == 'total':
            F1.append(indicator[i])

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
        indicator = list()
        indicator.append(precision)
        indicator.append(recall)
        indicator.append(f_score)
        indicator.append(mcc)
        if type == 'total':
            F1_gt.append(indicator[i])

# 读入gtf文件，并为其创建字典
def	gtf_read(dictgtf,path):
    file = open(path,"r")
    for line in file:
        seq = line.strip('\n').split('\t')
        chr = seq[0][3:]
        start = seq[3]
        end = seq[4]
        anno_type = seq[2]
        infor = seq[8]
        if chr not in dictgtf:
            dictgtf[chr] = list()
        dictgtf[chr].append([start,end,infor,anno_type])

# 对于INS 直接判断插入位点是否在gtf文件中即可（二分法）
def add_anotation_in(var_list, breakpoint):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = set()
    anno_type = set()
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
                anno_type.add(var_list[i][3])
        if breakpoint - int(var_list[i][0]) > 3000000:
            break
    return info,anno_type

# 对于其他类型的SV 需要判断overlap
def add_anotation_overlap(var_list, start, end):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = set()
    anno_type = set()
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
                anno_type.add(var_list[i][3])
        if start - int(var_list[i][0]) > 3000000:
            break
    for i in range(left + 1, len(var_list), 1):
        if int(var_list[i][0]) < end:
            if var_list[i][2] not in info:
                info.add(var_list[i][2])
                anno_type.add(var_list[i][3])
        else:
            break
    return info,anno_type


def gtf_call_read(dictgtf,path):
    vcf_reader = vcf.Reader(filename=path)
    annotation = list()
    annotation_type = list()
    for record in vcf_reader:
        if record.INFO['SVTYPE'] =='INS':
            if record.CHROM in dictgtf:
                a, b = add_anotation_in(dictgtf[record.CHROM], record.POS)
                annotation.append(a)
                annotation_type.append(b)
            else:
                annotation.append(set())
                annotation_type.append(set())

        else:
            if record.CHROM in dictgtf:
                a, b = add_anotation_overlap(dictgtf[record.CHROM],record.POS,record.INFO['END'])
                annotation.append(a)
                annotation_type.append(b)
            else:
                annotation.append(set())
                annotation_type.append(set())
        n = len(record.samples)

    file_handle = open(sys.argv[4], mode='w')
    i = 0
    f = codecs.open(path, mode='r', encoding='utf-8')
    line = f.readline()
    INFO = 0
    while line:
        a = line.split()
        if a[0][0] == "#" and a[0][2] != "I" and a[0][2] != "F":
            file_handle.write(line)
        elif a[0][0] == "#" and a[0][2] == "I":
            file_handle.write(line)
        elif a[0][0] == "#" and a[0][2] == "F" and INFO == 0:
            file_handle.write('##INFO=<ID=ANNOTYPE,Number=.,Type=String,Description=' + '"The annotation type"' + '>' + '\n')
            # file_handle.write('##INFO=<ID=DETAILS,Number=2,Type=String,Description=' + '"The details of annotation"' + '>' + '\n')
            file_handle.write(line)
            INFO = 1
        elif a[0][0] == "#" and a[0][2] == "F" and INFO != 0:
            file_handle.write(line)
        elif a[0][0] != "#":
            b = a[7] + ';'+'ANNOTYPE='

            if len(annotation_type[i])==0:
                b+='NULL'
            else:
                for x in annotation_type[i]:
                    b += x
                    b+=','
                b = b[:-1]
                b+=';'
            # b +=';'
            # b +='DETAILS='
            # for z in annotation[i]:
            #         b+=z

            i+=1

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

def print_txt(fx, svsize, F1):
    print("---------------------------------TOTAL---------------------------------")
    fx(total_perform, F1)
    print("---------------------------------INS---------------------------------")
    fx(svsize['INS'], F1)
    print("---------------------------------DEL---------------------------------")
    fx(svsize['DEL'], F1)
    print("---------------------------------DUP---------------------------------")
    fx(svsize['DUP'], F1)
    print("---------------------------------INV---------------------------------")
    fx(svsize['INV'], F1)
    print("---------------------------------BND---------------------------------")
    fx(svsize['BND'], F1)

def print_txt_convert(fx1,fx2,svsize,F1):
    print("---------------------------------TOTAL---------------------------------")
    fx1(total_perform, F1)
    print("---------------------------------INS---------------------------------")
    fx2('INS', F1)
    print("---------------------------------DEL---------------------------------")
    fx1(svsize['DEL'], F1)
    print("---------------------------------DUP---------------------------------")
    fx2('DUP', F1)
    print("---------------------------------INV---------------------------------")
    fx1(svsize['INV'], F1)
    print("---------------------------------BND---------------------------------")
    fx1(svsize['BND'], F1)

def photo_plot(dimension,F1,F1_gt):
    bar_width = 0.3
    index_presence = np.arange(len(dimension))
    index_gt = index_presence + bar_width
    plt.bar(index_presence, height=F1, width=bar_width, color='steelblue', label='presence')
    plt.bar(index_gt, height=F1_gt, width=bar_width, color='palevioletred', label='genotype')
    plt.legend()
    index = np.arange(len(dimension));
    for a, b in zip(index, F1):  # 柱子上的数字显示
        plt.text(a, b, '%.2f' % b, ha='center', va='bottom', fontsize=7);
    for a, b in zip(index + bar_width, F1_gt):
        plt.text(a, b, '%.2f' % b, ha='center', va='bottom', fontsize=7);
    #
    plt.xticks(index_presence + bar_width / 2, dimension)
    plt.ylabel('Indicator')
    plt.title('Indicator under different SV dimensions')
    plt.savefig('F1.png')

def comp_read(dictref,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        sv_type = record.INFO['SVTYPE']
        chr = record.CHROM
        try:
            genotype = record.samples[0]['GT']
        except:
            continue
        if sv_type not in dictref:
            dictref[sv_type] = dict()
            dictref[sv_type][chr]  = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        if sv_type == 'INS':
            dictref[sv_type][chr].append([record.POS,record.INFO['SVLEN'],genotype])
        elif sv_type == 'DEL' or sv_type == 'INV' or sv_type == 'DUP':
            dictref[sv_type][chr].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1, genotype])
        elif sv_type == "BND":
            form, chr, start, chr2, pos2 = phase_bnd(str(record.ALT[0]), chr, record.POS)
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
            dictref[sv_type][chr].append([start, chr2, pos2, form ,genotype])

def af_vcf_ref_read(dictref,path):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        # print(type(record.FILTER))
        sv_type = record.INFO['SVTYPE']
        chr = record.CHROM
        if len(record.FILTER) == 0:
            filt = 'PASS'
        else:
            filt = ''.join(str(v) for v in record.FILTER)
        try:
            genotype = record.samples[0]['GT']
        except:
            continue
        if sv_type not in dictref:
            dictref[sv_type] = dict()
            dictref[sv_type][chr]  = list()
        else:
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
        if sv_type == 'INS':
            dictref[sv_type][chr].append([record.POS,record.INFO['SVLEN'], record.CHROM,record.CHROM,record.CHROM, record.POS, record.ID, record.REF, ''.join(str(v) for v in record.ALT), record.QUAL, filt, record.INFO, record.FORMAT,record.samples, genotype])
        elif sv_type == 'DEL' or sv_type == 'INV' or sv_type == 'DUP':
            dictref[sv_type][chr].append([record.POS,record.INFO['END'],abs(record.INFO['END']-record.POS)+1, record.CHROM,record.CHROM, record.POS, record.ID, record.REF, ''.join(str(v) for v in record.ALT), record.QUAL, filt, record.INFO, record.FORMAT, record.samples,genotype])
        elif sv_type == "BND":
            form, chr, start, chr2, pos2 = phase_bnd(str(record.ALT[0]), chr, record.POS)
            if chr not in dictref[sv_type]:
                dictref[sv_type][chr] = list()
            dictref[sv_type][chr].append([start, chr2, pos2, form , record.CHROM, record.POS, record.ID, record.REF, ''.join(str(v) for v in record.ALT), record.QUAL, filt, record.INFO, record.FORMAT,record.samples,genotype])

def af_benchmark(dictref,dictcall,n):
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
                                if abs(call[0] - ref[0]) < threadshold[svtype]:
                                    if float(min(call[1],ref[1]) / max(call[1], ref[1])) >= bias:
                                        if ref[-n] == './.':
                                            ref[-n] = call[-1]
                                        break
                            elif svtype =='DEL' or svtype == 'INV' or svtype == 'DUP':
                                if max(ref[0] - threadshold[svtype], call[0]) < min(ref[1] + threadshold[svtype],call[1]):
                                    if float(min(call[2], ref[2]) /
                                     max(call[2], ref[2])) >= bias:
                                        if ref[-n] == './.':
                                            ref[-n] = call[-1]
                                        break
                            elif svtype == 'BND' :
                                if call[1] == ref[1] and abs(call[0] - ref[0]) < threadshold[svtype] and abs(call[2] - ref[2]) < threadshold[svtype] and call[3] == ref[3] :
                                    if ref[-n] == './.':
                                        ref[-n] = call[-1]
                                    break


if __name__ == '__main__':
    if sys.argv[1] == 'bench':
        if sys.argv[2] == '-h':
            print("\n")
            print("--------------------------------------usage---------------------------------------")
            print("usage: python main.py bench [base] [call] normal/convert [indicator] type/size")
            print("----------------------------------------------------------------------------------")
            print("\n")
        else:
            # print("输入ans的路径")
            refpath = sys.argv[2]
            dictref = dict()
            vcf_ref_read(dictref, refpath)

            # threadshold = 100  # 阈值
            bias = 0.7  # 长度偏移量

            # print("输入call的路径")
            callpath = sys.argv[3]
            dictcall = dict()
            vcf_ref_read(dictcall, callpath)

            svsize1 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen(),
                       'BND': init_stalen()}

            flag = 0
            benchmark(dictref, dictcall)
            statistics(dictref, dictcall, svsize1, -4)

            total_perform = init_stalen()
            for key in total_perform:
                for i in range(0, 6):
                    total_perform[key][i] = svsize1['INS'][key][i] + svsize1['DEL'][key][i] + \
                                            svsize1['DUP'][key][i] + svsize1['INV'][key][i] + svsize1['BND'][key][i]

            if sys.argv[4] == 'normal':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV', 'BND')
                dimension2 = ('total', '0-99', '100-499', '500-1k', '1k-5k', '5k-10k', '>10k')
                F1 = list()
                F1_gt = list()
                size_fig = list()
                size_fig_gt = list()
                print("No conversion between insertions and duplications")
                if sys.argv[5] == 'precision':
                    i = 0
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[5] == 'recall':
                    i = 1
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[5] == 'f_score':
                    i = 2
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[5] == 'mcc':
                    i = 3
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[6] == 'type':
                    photo_plot(dimension, F1, F1_gt)
                if sys.argv[6] == 'size':
                    x1 = size_fig[:7]
                    x2 = size_fig_gt[:7]
                    photo_plot(dimension2, x1, x2)

            if sys.argv[4] == 'convert':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV', 'BND')
                dimension2 = ('total', '0-99', '100-499', '500-1k', '1k-5k', '5k-10k', '>10k')
                F1 = list()
                F1_gt = list()
                size_fig = list()
                size_fig_gt = list()
                print("Converting between insertions and duplications")
                svsize2 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen(),
                           'BND': init_stalen()}
                ref_ins_call_dup()
                ref_dup_call_ins()
                statistics(dictref, dictcall, svsize2, -4)
                total_perform = init_stalen()
                for key in total_perform:
                    for i in range(0, 6):
                        total_perform[key][i] = svsize2['INS'][key][i] + svsize2['DEL'][key][i] + \
                                                svsize2['DUP'][key][i] + svsize2['INV'][key][i] + svsize1['BND'][key][i]
                for key in svsize2:
                    for type in svsize2[key]:
                        for i in range(0, 6):
                            svsize2[key][type][i] = svsize2[key][type][i] - svsize1[key][type][i]

                if sys.argv[5] == 'precision':
                    i = 0
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[5] == 'recall':
                    i = 1
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[5] == 'f_score':
                    i = 2
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[5] == 'mcc':
                    i = 3
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)

                if sys.argv[6] == 'type':
                    photo_plot(dimension, F1, F1_gt)
                if sys.argv[6] == 'size':
                    x1 = size_fig[:7]
                    x2 = size_fig_gt[:7]
                    photo_plot(dimension2, x1, x2)

    if sys.argv[1] == 'annotation':
        # print("输入gtf的路径")
        if sys.argv[2] == '-h':
            print("\n")
            print("--------------------------------------------------------usage--------------------------------------------------------")
            print("usage: python main.py annotation [gtf] [base] [output] [call] [anno_type] normal/convert [indicator] size/type")
            print("----------------------------------------------------------------------------------------------------------------------")
            print("\n")
        else:
            gtfpath = sys.argv[2]
            dictgtf = dict()
            gtf_read(dictgtf, gtfpath)
            print("annotating")
            gtfbasepath = sys.argv[3]
            gtf_call_read( dictgtf, gtfbasepath)

            path = sys.argv[4]
            dictref = dict()
            gtf_vcf_ref_read(dictref, path)
            # print(dicttest)

            call_path = sys.argv[5]
            dictcall = dict()
            vcf_ref_read(dictcall, call_path)

            # threadshold = 100  # 阈值
            bias = 0.7  # 长度偏移量
            svsize1 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen(),
                       'BND': init_stalen()}

            flag = 1
            anno_type = sys.argv[6]
            benchmark(dictref, dictcall)
            statistics(dictref, dictcall, svsize1, -5)

            total_perform = init_stalen()
            for key in total_perform:
                for i in range(0, 6):
                    total_perform[key][i] = svsize1['INS'][key][i] + svsize1['DEL'][key][i] + \
                                            svsize1['DUP'][key][i] + svsize1['INV'][key][i] + svsize1['BND'][key][i]

            if sys.argv[7] == 'normal':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV', 'BND')
                dimension2 = ('total', '0-99', '100-499', '500-1k', '1k-5k', '5k-10k', '>10k')
                F1 = list()
                F1_gt = list()
                size_fig = list()
                size_fig_gt = list()

                if sys.argv[8] == 'precision':
                    i = 0
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[8] == 'recall':
                    i = 1
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[8] == 'f_score':
                    i = 2
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)
                if sys.argv[8] == 'mcc':
                    i = 3
                    print_txt(print_result, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt(print_result_gt, svsize1, F1_gt)

                if sys.argv[9] == 'type':
                    photo_plot(dimension, F1, F1_gt)
                if sys.argv[9] == 'size':
                    x1 = size_fig[:7]
                    x2 = size_fig_gt[:7]
                    photo_plot(dimension2, x1, x2)

            if sys.argv[7] == 'convert':
                dimension = ('TOTAL', 'INS', 'DEL', 'DUP', 'INV', 'BND')
                dimension2 = ('total', '0-99', '100-499', '500-1k', '1k-5k', '5k-10k', '>10k')
                F1 = list()
                F1_gt = list()
                size_fig = list()
                size_fig_gt = list()
                print("Converting between insertions and duplications")
                svsize2 = {'INS': init_stalen(), 'DEL': init_stalen(), 'DUP': init_stalen(), 'INV': init_stalen(),
                           'BND': init_stalen()}
                ref_ins_call_dup()
                ref_dup_call_ins()
                statistics(dictref, dictcall, svsize2, -5)
                total_perform = init_stalen()
                for key in total_perform:
                    for i in range(0, 6):
                        total_perform[key][i] = svsize2['INS'][key][i] + svsize2['DEL'][key][i] + \
                                                svsize2['DUP'][key][i] + svsize2['INV'][key][i] + svsize1['BND'][key][i]
                for key in svsize2:
                    for type in svsize2[key]:
                        for i in range(0, 6):
                            svsize2[key][type][i] = svsize2[key][type][i] - svsize1[key][type][i]

                if sys.argv[8] == 'precision':
                    i = 0
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[8] == 'recall':
                    i = 1
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[8] == 'f_score':
                    i = 2
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)
                if sys.argv[8] == 'mcc':
                    i = 3
                    print_txt_convert(print_result, print_result_addition, svsize1, F1)
                    print("\n")
                    print("-----------------------------------GT-----------------------------------")
                    print_txt_convert(print_result_gt, print_result_addition_gt, svsize1, F1_gt)


                if sys.argv[9] == 'type':
                    photo_plot(dimension, F1, F1_gt)
                if sys.argv[9] == 'size':
                    x1 = size_fig[:7]
                    x2 = size_fig_gt[:7]
                    photo_plot(dimension2, x1, x2)

    if sys.argv[1] == 'af':
        if sys.argv[2] == '-h':
            print("\n")
            print("-------------------------------usage-------------------------------")
            print("usage: python main.py af [base] [call1] [call2] ... [calln]")
            print("-------------------------------------------------------------------")
            print("\n")
        else:
            dictref = dict()
            ref_path = sys.argv[2]
            af_vcf_ref_read(dictref, ref_path)

            # threadshold = 100  # 阈值
            bias = 0.7  # 长度偏移量

            for svtype in dictref:
                for chr in dictref[svtype]:
                    for sv in dictref[svtype][chr]:
                        for i in range(0, len(sys.argv) - 3):
                            sv.append('./.')

            for i in range(0, len(sys.argv) - 3):
                print("compare with call"+str(i+1))
                dictcall = dict()
                call_path = sys.argv[i + 3]
                comp_read(dictcall, call_path)
                # print(1 + i)
                af_benchmark(dictref, dictcall, 1 + i)

            n = len(sys.argv) - 2
            # print(n)

            for svtype in dictref:
                for chr in dictref[svtype]:
                    for sv in dictref[svtype][chr]:
                        AF = 0
                        for i in range(0, len(sys.argv) - 2):
                            if sv[-(i + 1)] == '1/1':
                                AF += 2
                            if sv[-(i + 1)] == '0/1':
                                AF += 1
                        AF = AF / (2 * n)
                        sv.append(AF)

            head_file = open(sys.argv[2], 'r')
            line = head_file.readline()
            file = open('gene_af.vcf', 'w')
            INFO = 0
            while line:
                a = line.split()
                if a[0][0] == "#" and a[0][2] != "I" and a[0][2] != "F":
                    if a[0][1] == "C":
                        file.write(line.strip('\n'))
                        file.write('\t')
                        for i in range(0, len(sys.argv) - 3):
                            file.write("Sample" + str(i + 1) + '\t')
                        file.write('\n')
                    else:
                        file.write(line)

                elif a[0][0] == "#" and a[0][2] == "I":
                    file.write(line)
                elif a[0][0] == "#" and a[0][2] == "F" and INFO == 0:
                    file.write('##INFO=<ID=AF,Number=1,Type=Float,Description=' + '"Allele frequency"' + '>' + '\n')
                    file.write(line)
                    INFO = 1
                elif a[0][0] == "#" and a[0][2] == "F" and INFO != 0:
                    file.write(line)
                line = head_file.readline()

            for svtype in dictref:
                for chr in dictref[svtype]:
                    for sv in dictref[svtype][chr]:
                        for i in range(4, 11):
                            file.write(str(sv[i]) + '\t')
                        for key in sv[11]:
                            file.write(key)
                            if isinstance(sv[11][key], list):
                                file.write('=')
                                file.write(','.join(str(v) for v in sv[11][key]))
                                file.write(';')
                            if isinstance(sv[11][key], bool):
                                file.write(';')
                            elif not isinstance(sv[11][key], list):
                                file.write('=')
                                file.write(str(sv[11][key]))
                                file.write(';')
                        file.write('AF=' + str(round(sv[-1], 4)))
                        file.write('\t')
                        file.write('GT' + '\t')
                        for i in range(0, len(sys.argv) - 2):
                            file.write(str(sv[-(len(sys.argv) - i - 1)]) + '\t')
                        file.write('\n')