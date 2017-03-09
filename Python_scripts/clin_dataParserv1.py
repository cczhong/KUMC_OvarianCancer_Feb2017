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
subtypeData.outputProlifKeys()


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
#print('row#: ', len(expressionMatrix), '\n col#: ', len(expressionMatrix[0]))

for i in range(0, len(expressionMatrix[1])):
    if platStatusDict.get(expressionMatrix[0][i]) != 'None':
        if fallopian_Dict.get(expressionMatrix[0][i]) == 'Sensitive' or fallopian_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i])
            expVector.append(fallopian_Dict.get(expressionMatrix[0][i]))
            fallopian_exp_data.append(expVector)
        elif proliferative_Dict.get(expressionMatrix[0][i]) == 'Sensitive' or proliferative_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i])
            expVector.append(proliferative_Dict.get(expressionMatrix[0][i]))
            prolif_exp_data.append(expVector)
        elif mesenchymal_Dict.get(expressionMatrix[0][i]) == 'Sensitive' or mesenchymal_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i])
            expVector.append(mesenchymal_Dict.get(expressionMatrix[0][i]))
            mesen_exp_data.append(expVector)
        elif immunoreactive_Dict.get(expressionMatrix[0][i]) == 'Sensitive' or immunoreactive_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            expVector = []
            expVector.append(expressionMatrix[0][i])
            expVector.append(immunoreactive_Dict.get(expressionMatrix[0][i])) 
            immuno_exp_data.append(expVector)

platStat = ['Platinum Status']
final_out_matrix = []

for i in range(0, len(expressionMatrix)):
    final_out_matrix.append(['']*len(expressionMatrix[1]))

#print('Length of final output matrix: ', len(final_out_matrix), '\nNumber of rows in matrix: ', len(final_out_matrix[0]))

for i in range(1, len(expressionMatrix)):
    platStat.append([]*len(expressionMatrix[0]))
    


phenotypeMatrix = []

for i in range(0, len(expressionMatrix[0])):
    phenotypeMatrix.append([' ']*3)
#print(len(phenotypeMatrix), phenotypeMatrix[0])
subtypeCol = ['Subtype']


for i in range(0, len(expressionMatrix[0])):
    if platStatusDict.get(expressionMatrix[0][i]) != 'None':
        if fallopian_Dict.get(expressionMatrix[0][i]) == 'Sensitive':
            phenotypeMatrix[i][0] = expressionMatrix[0][i]
            phenotypeMatrix[i][1] = 'Fallopian'
            phenotypeMatrix[i][2] = 'Sensitive'
        elif fallopian_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            phenotypeMatrix[i][0] = expressionMatrix[0][i]
            phenotypeMatrix[i][1] = 'Fallopian' 
            phenotypeMatrix[i][2] = 'Resistant'
        elif proliferative_Dict.get(expressionMatrix[0][i]) == 'Sensitive':
            phenotypeMatrix[i][0] = expressionMatrix[0][i]
            phenotypeMatrix[i][1] = 'Proliferative'
            phenotypeMatrix[i][2] = 'Sensitive'
            platStat[i] = 'Sensitive'
        elif proliferative_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            phenotypeMatrix[i][0] = expressionMatrix[0][i] 
            phenotypeMatrix[i][1] ='Proliferative'
            phenotypeMatrix[i][2] = 'Resistant'
            platStat[i] = 'Resistant'
        elif mesenchymal_Dict.get(expressionMatrix[0][i]) == 'Sensitive':
            phenotypeMatrix[i][0] = expressionMatrix[0][i] 
            phenotypeMatrix[i][1] = 'Mesenchymal'
            phenotypeMatrix[i][2] = 'Sensitive'
        elif mesenchymal_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            phenotypeMatrix[i][0] = expressionMatrix[0][i] 
            phenotypeMatrix[i][1] = 'Mesenchymal'
            phenotypeMatrix[i][2] = 'Resistant'
        elif immunoreactive_Dict.get(expressionMatrix[0][i]) == 'Sensitive':
            phenotypeMatrix[i][0] = expressionMatrix[0][i] 
            phenotypeMatrix[i][1] = 'Immunoreactive'
            phenotypeMatrix[i][2] = 'Sensitive'
        elif immunoreactive_Dict.get(expressionMatrix[0][i]) == 'Resistant':
            phenotypeMatrix[i][0] = expressionMatrix[0][i] 
            phenotypeMatrix[i][1] = 'Immunoreactive'
            phenotypeMatrix[i][2] = 'Resistant'
        else:
            if i == 0:
                phenotypeMatrix[i][0] = 'Sample_ID '
                phenotypeMatrix[i][1] = 'Subtype'
                phenotypeMatrix[i][2] = 'Platinum Status'
            else:
                phenotypeMatrix[i][0] = expressionMatrix[0][i]
                phenotypeMatrix[i][1] = 'Unknown'
                phenotypeMatrix[i][2] = 'Unknown'

phenotypeOut = open('phenotype_matrix.txt', 'w+')
for i in range(0, len(phenotypeMatrix)):
    for j in range(0, len(phenotypeMatrix[0])):
        phenotypeOut.write(phenotypeMatrix[i][j])
        phenotypeOut.write('\t')
    phenotypeOut.write('\n')

phenotypeOut.close()






