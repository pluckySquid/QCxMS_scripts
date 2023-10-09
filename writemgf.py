from pyteomics.mgf import write
import csv



file = open("result.csv", "r")
data = list(csv.reader(file, delimiter=","))
mass = [row[0] for row in data]
mass = [float(x) for x in mass]
ratio = [row[1] for row in data]
ratio = [float(x) for x in ratio]
print(mass)

ans_mass = []
ans_ratio = []

for i in range(0, len(ratio)):
	if ratio[i] > 0.1:
		ans_mass.append(mass[i])
		ans_ratio.append(ratio[i])
#write([{'m/z array':[1,2,3,4],'intensity array':[0.1, 0.2, 0.4, 1],'params':{'scan':1}}], output_path)
write([{'m/z array':mass,'intensity array':ratio,'params':{'scan':1}}], "result.mgf")
write([{'m/z array':ans_mass,'intensity array':ans_ratio,'params':{'scan':1}}], "ans_result.mgf")
