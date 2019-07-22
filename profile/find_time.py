import numpy as np
import matplotlib.pyplot as plt

filename1 = "/Users/levius/Desktop/research_genomic_variant_call/my_method/Variant_Calling/profile/baseline_profile"
filename2 = "/Users/levius/Desktop/research_genomic_variant_call/my_method/Variant_Calling/profile/prune_profile"

filter_base_time = 0
genotype_base_time = 0
region_base_num = 0

filter_prune_time = 0
genotype_prune_time = 0
region_prune_time = 0

original_work = []
filter_work = []
redo_work = []

file1 = open(filename1, 'r')
line = file1.readline()
line = file1.readline()
while line != '':
    region_base_num += 1
    line = file1.readline()
    line = line.split()
    filter_base_time += int(line[-1])
    line = file1.readline()
    line = line.split()
    genotype_base_time += int(line[-1])
    line = file1.readline()

file1.close()
print("#base region:", region_base_num)
print("filter base time(ms):", filter_base_time)
print("genotype base time(ms):", genotype_base_time)
print("baseline time after pairHMM(ms):", genotype_base_time + filter_base_time)

file2 = open(filename2, 'r')
line = file2.readline()
line = file2.readline()
while line != '':
    region_prune_time += 1
    line = file2.readline()
    line = line.split()
    genotype_prune_time += int(line[-1])
    line = file2.readline()
    line = line.split()
    original_work.append(int(line[-1]))
    line = file2.readline()
    line = line.split()
    filter_work.append(int(line[-1]))
    line = file2.readline()
    line = line.split()
    redo_work.append(int(line[-1]))
    line = file2.readline()

file2.close()
print("#prune region:", region_prune_time)
# print("filter prune time(ms):", filter_prune_time)
print("genotype prune time(ms):", genotype_prune_time)
print("prune time after pairHMM(ms):", genotype_prune_time + filter_prune_time)

filter_work = np.array(filter_work)
original_work = np.array(original_work)
redo_work = np.array(redo_work)

total_work = original_work.sum()
total_filter = filter_work.sum()
total_redo = redo_work.sum()

print("total original work(cell):", total_work)
print("total filter work(cell):", total_filter)
print("total redo work(cell):", total_redo)

print("filter percentage:", total_filter/total_work)
print("redo percentage:", total_redo/total_work)