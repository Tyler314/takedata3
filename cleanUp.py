from sys import argv

f = open(argv[1], 'r')
if len(argv) > 2:
    newFile = argv[2]
else:
    newFile = argv[1].split('.')[0] + '_clean_up.' + argv[1].split('.')[1]

g = open(newFile, 'w')

for line in f:
    g.write(line.replace('\t', 8*' ').rstrip() + '\n')

f.close()
g.close()
