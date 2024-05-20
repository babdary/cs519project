import random
import sys

nodes=int(sys.argv[1])
edges_per_node=int(sys.argv[2])
f=open('inputMatrix', 'w')
f.write(str(nodes) + '\n')
f.write(str(edges_per_node*nodes)+ '\n')

for i in range(nodes):
    for j in range(edges_per_node):
        f.write(str(int((i+1)*random.random()+1)) +' ' + str(int((j+1)*random.random()+1))+ ' '+ str(int(100*random.random()))+ '\n')

f.close()