import csv
import sys
from subtypeIdClasses import subtypeIdClasses

def printDataMatrix(matrix):
    for i in range (len(matrix)):
        for j in range(0, 18):
            sys.stdout.write(matrix[i][j])
            sys.stdout.write('\t')
        sys.stdout.write('\n')

def prelimPlatStatParse(matrix, s_type, r_type, early, blank):
    for i in range(len(matrix)):
        if matrix[i][1] == 'Sensitive':
            s_type.append(matrix[i])
        elif matrix[i][1] == 'Resistant':
            r_type.append(matrix[i])
        elif matrix[i][1] == 'Tooearly':
            early.append(matrix[i])
        elif matrix[i][1] == '':
            blank.append(matrix[i])
            

def fullPlatStatParse(matrix, s_type, r_type):
    for i in range(len(matrix)):
        if matrix[i][1] == 'Sensitive':
            s_type.append(matrix[i])
        elif matrix[i][1] == 'Resistant':
            r_type.append(matrix[i])
        elif matrix[i][1] == 'Tooearly':
            if matrix[i][2] == 'COMPLETE RESPONSE':
                s_type.append(matrix[i])
            elif matrix[i][2] == 'PARTIAL RESPONSE':
                s_type.append(matrix[i])
            elif matrix[i][2] == 'STABLE DISEASE':
                r_type.append(matrix[i])
            elif matrix[i][2] == 'PROGRESSIVE DISEASE':
                r_type.append(matrix[i])
        elif matrix[i][1] == '':
            if matrix[i][2] == 'COMPLETE RESPONSE':
                s_type.append(matrix[i])
            elif matrix[i][2] == 'PARTIAL RESPONSE':
                s_type.append(matrix[i])
            elif matrix[i][2] == 'STABLE DISEASE':
                r_type.append(matrix[i])
            elif matrix[i][2] == 'PROGRESSIVE DISEASE':
                r_type.append(matrix[i])

def countResistAndPartialResponse(matrix):
    count = 0
    for i in range(len(matrix)):
        if matrix[i][1] == 'Resistant' and matrix[i][2] == 'PARTIAL RESPONSE':
            count +=1
    return count

def countResistComplete(matrix):
    count = 0
    for i in range(len(matrix)):
        if matrix[i][1] == 'Resistant' and matrix[i][2] == 'COMPLETE RESPONSE':
            count += 1
    return count

def searchMatrices(searchWord, m1, m2, m3, m4, m5, m6, m7, m8):
    for i in range(0, len(m1)):
        pass
dataList = list(csv.reader(open('ov_tcga_pub_clinical_data.tsv'), delimiter='\t'))


#will contain all samples of patient data in tuples of the form (sample ID, platinum status, patient response)
fallopian_st_matrix = [] 
mesenchymal_st_matrix = []
proliferative_st_matrix = []
immunoreactive_st_matrix = []

r_fallopian_st = []
r_proliferative_st = []
r_mesenchymal_st = []
r_immunoreactive_st = []


s_fallopian_st = []
s_proliferative_st = []
s_mesenchymal_st = []
s_immunoreactive_st = []

e_fallopian_st = []
e_proliferative_st = []
e_mesenchymal_st = []
e_immunoreactive_st = []

b_fallopian_st = []
b_proliferative_st = []
b_mesenchymal_st = []
b_immunoreactive_st = []

subtypeData = subtypeIdClasses()


#groups data table by subtype using sample IDs.
for i in range(len(dataList)):
    if subtypeData.bool_searchFallopianIds_d(dataList[i][2]) == True:
        addId = subtypeData.get_fallopianIndex(dataList[i][2])
        temp =[subtypeData.p_fallopianIdKeys[addId], dataList[i][0], dataList[i][1]]
        fallopian_st_matrix.append(temp)
        
    elif subtypeData.bool_searchProlifIds_d(dataList[i][2]) == True:
        addId = subtypeData.get_prolifIndex(dataList[i][2])
        temp =[subtypeData.p_proliferativeIdKeys[addId], dataList[i][0], dataList[i][1]]
        proliferative_st_matrix.append(temp)
        
    elif subtypeData.bool_searchMesenIds_d(dataList[i][2]) == True:
        addId = subtypeData.get_mesenIndex(dataList[i][2])
        temp =[subtypeData.p_mesenchymalIdKeys[addId], dataList[i][0], dataList[i][1]]
        mesenchymal_st_matrix.append(temp)
        
    elif subtypeData.bool_searchImmunoIds_d(dataList[i][2]) == True:
        addId = subtypeData.get_immunoIndex(dataList[i][2])
        temp =[subtypeData.p_immunoreactiveIdKeys[addId], dataList[i][0], dataList[i][1]]
        immunoreactive_st_matrix.append(temp)



fullPlatStatParse(fallopian_st_matrix, s_fallopian_st, r_fallopian_st)
fullPlatStatParse(proliferative_st_matrix, s_proliferative_st, r_proliferative_st)
fullPlatStatParse(mesenchymal_st_matrix, s_mesenchymal_st, r_mesenchymal_st)
fullPlatStatParse(immunoreactive_st_matrix, s_immunoreactive_st, r_immunoreactive_st)
platStatusDict = {}
fallopian_Dict = {}
proliferative_Dict = {}
mesenchymal_Dict = {}
immunoreactive_Dict = {}

for i in range(0, len(s_fallopian_st)):
    platStatusDict[s_fallopian_st[i][0]] = 'Sensitive'
    fallopian_Dict[s_fallopian_st[i][0]] = 'Sensitive'
for i in range(0, len(r_fallopian_st)):
    platStatusDict[r_fallopian_st[i][0]] = 'Resistant'
    fallopian_Dict[r_fallopian_st[i][0]] = 'Resistant'
for i in range(0, len(s_proliferative_st)):
    platStatusDict[s_proliferative_st[i][0]] = 'Sensitive'
    proliferative_Dict[s_proliferative_st[i][0]] = 'Sensitive'
    
for i in range(0, len(r_proliferative_st)):
    platStatusDict[r_proliferative_st[i][0]] = 'Resistant'
    proliferative_Dict[r_proliferative_st[i][0]] = 'Resistant'
for i in range(0, len(s_mesenchymal_st)):
    platStatusDict[s_mesenchymal_st[i][0]] = 'Sensitive'
    mesenchymal_Dict[s_mesenchymal_st[i][0]] = 'Sensitive'
for i in range(0, len(r_mesenchymal_st)):
    platStatusDict[r_mesenchymal_st[i][0]] = 'Resistant'
    mesenchymal_Dict[r_mesenchymal_st[i][0]] = 'Resistant'
for i in range(0, len(s_immunoreactive_st)):
    platStatusDict[s_immunoreactive_st[i][0]] = 'Sensitive'
    immunoreactive_Dict[s_immunoreactive_st[i][0]] = 'Sensitive'
for i in range(0, len(r_immunoreactive_st)):
    platStatusDict[r_immunoreactive_st[i][0]] = 'Resistant' 
    immunoreactive_Dict[r_immunoreactive_st[i][0]] = 'Resistant'




fallopian_exp_data = []
prolif_exp_data = []
mesen_exp_data = []
immuno_exp_data = []

#Parsing microarray data
expressionMatrix = list(csv.reader(open('TCGA_489_UE.txt'), delimiter = '\t'))


for i in range(1, len(expressionMatrix[1])):
    if platStatusDict.get(expressionMatrix[0][i-1]) != 'None':
        if fallopian_Dict.get(expressionMatrix[0][i-1]) == 'Sensitive' or fallopian_Dict.get(expressionMatrix[0][i-1]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i-1])
            expVector.append(fallopian_Dict.get(expressionMatrix[0][i-1]))
            padStr = ''
            padVector = []
            for rowInd in range(1, 11865):
                expVector.append(expressionMatrix[rowInd][i])
                padStr = 'Column ' + str(i)
                padVector.append(padStr)
            if i == 1:
                fallopian_exp_data.append(padVector)
            fallopian_exp_data.append(expVector)
        elif proliferative_Dict.get(expressionMatrix[0][i-1]) == 'Sensitive' or proliferative_Dict.get(expressionMatrix[0][i-1]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i-1])
            expVector.append(proliferative_Dict.get(expressionMatrix[0][i-1]))
            padStr = ''
            padVector = []
            for rowInd in range(1, 11865):
                expVector.append(expressionMatrix[rowInd][i])
                padStr = 'Column ' + str(i)
                padVector.append(padStr)
            if i == 1:
                prolif_exp_data.append(padVector)
            prolif_exp_data.append(expVector)
        elif mesenchymal_Dict.get(expressionMatrix[0][i-1]) == 'Sensitive' or mesenchymal_Dict.get(expressionMatrix[0][i-1]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i-1])
            expVector.append(mesenchymal_Dict.get(expressionMatrix[0][i-1]))
            padStr = ''
            padVector = []
            for rowInd in range(1, 11865):
                padStr= 'Column ' + str(i)
                padVector.append(padStr)
                expVector.append(expressionMatrix[rowInd][i])
            if i == 1:
                mesen_exp_data.append(padVector)
            mesen_exp_data.append(expVector)
        elif immunoreactive_Dict.get(expressionMatrix[0][i-1]) == 'Sensitive' or immunoreactive_Dict.get(expressionMatrix[0][i-1]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i-1])
            expVector.append(immunoreactive_Dict.get(expressionMatrix[0][i-1]))
            padStr = ''
            padVector = []
            for rowInd in range(1, 11865):
                expVector.append(expressionMatrix[rowInd][i])
                padStr = 'Column ' + str(i)
                padVector.append(padStr)
            if i == 1:
                immuno_exp_data.append(padVector) 
            immuno_exp_data.append(expVector)



#data output
f_out = open('fallopian.txt', 'w+')
p_out = open('proliferative.txt', 'w+')
m_out = open('mesenchymal.txt', 'w+')
i_out = open('immunoreactive.txt','w+')

#write1 = csv.writer(f_out,  quoting=csv.QUOTE_ALL)
#write2 = csv.writer(p_out,  quoting=csv.QUOTE_ALL)
#write3 = csv.writer(m_out,  quoting=csv.QUOTE_ALL)
#write4 = csv.writer(i_out, quoting=csv.QUOTE_ALL)

for i in range(0, len(fallopian_exp_data)):
    for j in range(0, len(fallopian_exp_data[0])):
        f_out.write(fallopian_exp_data[i][j])
        f_out.write('\t')
    f_out.write('\n')
for i in range(0, len(prolif_exp_data)):
    for j in range(0, len(prolif_exp_data[0])):
        p_out.write(prolif_exp_data[i][j])
        p_out.write('\t')
    p_out.write('\n')
for i in range(0, len(mesen_exp_data)):
    for j in range(0, len(mesen_exp_data[0])):
        m_out.write(mesen_exp_data[i][j])
        m_out.write('\t')
    m_out.write('\n')
for i in range(0, len(immuno_exp_data)):
    for j in range(0, len(immuno_exp_data[0])):
        i_out.write(immuno_exp_data[i][j])
        i_out.write('\t')
    i_out.write('\n')


f_out.close()
p_out.close()
m_out.close()
i_out.close()




