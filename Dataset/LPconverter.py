from os import listdir
from os.path import isfile, join


#Reading the TXT file

#f = open("scp41.txt", "r")
mypath="./AA/"#"./Rail/" #"./Original/"
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files.sort()
print(files)
for k in files:
    rail = False
    AA = True
    if(not rail and not AA):
        with open(mypath+k) as f:
            flat_list=[int(word) for line in f for word in line.split()]
        #line = f.readline().split(' ')
        #print(line)
        #rows=line[1]
        #columns=line[2]

        #print(flat_list)
        # Reading the number of rows and columns
        rows = flat_list[0]
        columns = flat_list[1]
        print("rows ",rows, " columns ", columns)
    
        #Reading the costs of the variables(columns)
        costs=flat_list[2:(columns+2)]
        #print(costs, "\nCosts length ", len(costs))

        constraints = []
        index=columns+2

        for i in range(rows):
            constraint = flat_list[index+1:index+1+flat_list[index]] 
            constraints.append(constraint)
            #print("constraint ", i,", it has ", flat_list[index], " elements \n",constraint)
            index+=flat_list[index]+1
    if(rail):
        costs=[]
        f = open(mypath+k, 'r')
        counter=1 
        first = True
        for line in f:
            line = line.split()
            if first == True:
                rows = int(line[0])
                columns = int(line[1])
                constraints = [[] for _ in range(rows)]
                #print("rows ",rows, " columns ", columns)
                first=False
            else:
                #print("managing constraint ", counter)
                costs.append(int(line[0]))
                for i in range(int(line[1])):
                    #print(line[2+i])
                    constraints[int(line[2+i])-1].append(counter)
                counter+=1
    if(AA):
        with open(mypath+k) as f:
            flat_list=[int(word) for line in f for word in line.split()]
        
        rows = flat_list[1]
        columns = flat_list[0]
        print("rows ",rows, " columns ", columns)
        costs=flat_list[2:(columns+2)]
       
        constraints = [[] for _ in range(rows)]
        print("constraints = ", len(constraints))
        index=columns+2
        for i in range(columns):
            print("column ", i, " appears in ", flat_list[index], " constraints")
            for j in range(1, flat_list[index]+1):
                constraints[flat_list[index+j]-1].append(i+1)
                
            index+=flat_list[index]+1
        print(constraints)
    


    # writing the corresponding LP format
    new_name="./LPAA/"+k[0:-3]+"lp" #"./LPRail/"+k[0:-3]+"lp"
    f = open(new_name, "w")
    f.write("Minimize\n obj: ")
    linelength = 6 
    for i in range(columns-1):
        element = str(costs[i]) + " x" + str(i + 1)
        linelength+=len(element)
        if linelength + 3 < 80:
            f.write(element+" + ")
        else:
            f.write(element+"\n\t  + ")
            linelength = 6

    f.write(str(costs[columns-1]) + " x" + str(columns)+"\nSubject To:\n")

    linelength = 6 
    for i in range(rows):
        f.write(" c"+str(i+1)+":\t")
        for j in range(len(constraints[i])-1):
            element = "x" + str(constraints[i][j])
            linelength+=len(element)
            if linelength + 3 < 70:
                f.write(element+" + ")
            else:
                f.write(element+"\n\t  + ")
                linelength = 6
        print(i)
        f.write("x" + str(constraints[i][len(constraints[i])-1])+" >= 1\n")
        linelength = 6

    f.write("Bounds\n")
    for i in range(columns):
        f.write(" 0 <= x"+str(i+1)+ " <= 1\n")

    f.write("Binaries\n")

    linelength = 0
    for i in range(columns):
        element = " x"+str(i+1)
        linelength+=len(element)
        if linelength <= 80:
            f.write(element)
        else:
            f.write(element+"\n")
            linelength = len(element)

    f.write("\nEnd\n")
            



#print(f.read())
