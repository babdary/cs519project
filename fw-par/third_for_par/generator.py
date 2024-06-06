import random

random.seed(5)
size = 10

output = []

for x in range(size):
    temp = []
    for x in range(size):
        temp.append([])
    output.append(temp)

max_length = 1

for x in range(size):
    for y in range(x, size):
        if(x == y):
            #print(output)
            output[x][y] = 0
        else:
            number = int(random.random()*10)+1
            roll = int(random.random()*4)
            #print(roll)
            if(roll == 0):
                output[x][y] = 99999
                output[y][x] = 99999
            elif(roll == 1):
                output[x][y] = int(number)
                output[y][x] = 99999
            elif(roll == 2):
                output[x][y] = 99999
                output[y][x] = int(number)
            elif(roll == 3):
                output[x][y] = int(number)
                output[y][x] = int(number)

text1 = "{ "
for x in output:
    text1 += "{ "
    text1 += str(x)[1:-1]
    text1 += "}, "
text1 = text1[:-2]
text1 += "};"
#print(text1)
f = open("D://work//Classes//CS519//project//matrix.txt", "w")
f.write(text1)
f.close()

#output_linear = []
#for x in range(size):
#    for y in range(size):
#        output_linear.append(output[x][y])
#print("{ "+str(output_linear)[1:-1]+" };")
