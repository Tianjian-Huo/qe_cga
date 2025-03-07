"""
write by fuli 2022.10.10

before using this file, please do 'tail -28 dmol.outmol|head -16 > str.txt' first
This file is to extract structure from test.txt
find optimized structure and write it into .car file

"""

title = "!BIOSYM archive 3 \nPBC=OFF \ndmol opted structure  \n!DATE     2022-10-10\n"

filename = "str.txt"
new_file = "test.car"

N = 0

f_new = open(new_file, "w")
f_new.write(title)
f = open(filename)
while True:
    N += 1
    line = f.readline()
    if line == '':
        break
    XYZ = line.split()
    print(XYZ[-3],XYZ[-2],XYZ[-1])
    f_new.write("Ag" + str(N) + "\t" + XYZ[-3] + "\t" + XYZ[-2] + "\t" + XYZ[-1] + "   XXXX 1      xx      Ag  0.000\n")
f_new.write("end\nend\n")
f_new.close()
