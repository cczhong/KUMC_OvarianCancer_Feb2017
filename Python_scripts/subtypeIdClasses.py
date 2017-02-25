import sys
class subtypeIdClasses:
    def __init__(self):    
        self.d_fallopianIdKeys = []
        self.p_fallopianIdKeys = []
        
        self.d_proliferativeIdKeys = []
        self.p_proliferativeIdKeys = []
        
        self.d_mesenchymalIdKeys = []
        self.p_mesenchymalIdKeys = []
        
        self.d_immunoreactiveIdKeys = []
        self.p_immunoreactiveIdKeys = []
        
        self.initKeyVectors()
        
    def initKeyVectors(self):
        with open(sys.argv[1], 'r') as f:
            for line in f:
                readString = ''
                readString += line[0:20]
                self.d_fallopianIdKeys.append(readString)
                self.p_fallopianIdKeys.append(self.swapDashes(readString))

        with open(sys.argv[2], 'r') as p:
            for line in p:
                readP = ''
                readP += line[0:20]
                self.d_proliferativeIdKeys.append(readP)
                self.p_proliferativeIdKeys.append(self.swapDashes(readP))
        with open(sys.argv[3], 'r') as m:
            for line in m:
                readM = ''
                readM += line[1:21]
                self.d_mesenchymalIdKeys.append(readM)
                self.p_mesenchymalIdKeys.append(self.swapDashes(readM))
        with open(sys.argv[4], 'r') as im:
            for line in im:
                readI = ''
                readI += line[1:21]
                self.d_immunoreactiveIdKeys.append(readI)
                self.p_immunoreactiveIdKeys.append(self.swapDashes(readI))
        
        
    def bool_searchFallopianIds_d(self,compId):
        for i in range(len(self.d_fallopianIdKeys)):
            if self.d_fallopianIdKeys[i][:15] == compId:
                return True
    def bool_searchProlifIds_d(self, compId):
        for i in range(len(self.d_proliferativeIdKeys)):
            if self.d_proliferativeIdKeys[i][:15] == compId:
                return True
    def bool_searchImmunoIds_d(self, compId):
        for i in range(len(self.d_immunoreactiveIdKeys)):
            if self.d_immunoreactiveIdKeys[i][:15] == compId:
                return True
    def bool_searchMesenIds_d(self, compId):
        for i in range(len(self.d_mesenchymalIdKeys)):
            if self.d_mesenchymalIdKeys[i][:15] == compId:
                return True  
            
    def bool_searchFallopianIds_p(self,compId):
        for i in range(len(self.p_fallopianIdKeys)):
            if self.p_fallopianIdKeys[i][:15] == compId:
                return True
    def bool_searchProlifIds_p(self, compId):
        for i in range(len(self.p_proliferativeIdKeys)):
            if self.p_proliferativeIdKeys[i][:15] == compId:
                return True
    def bool_searchImmunoIds_p(self, compId):
        for i in range(len(self.p_immunoreactiveIdKeys)):
            if self.p_immunoreactiveIdKeys[i][:15] == compId:
                return True
    def bool_searchMesenIds_p(self, compId):
        for i in range(len(self.p_mesenchymalIdKeys)):
            if self.p_mesenchymalIdKeys[i][:15] == compId:
                return True              
                  
    def get_fallopianIndex(self, compId):
        for i in range(len(self.d_fallopianIdKeys)):
            if self.d_fallopianIdKeys[i][:15] == compId:
                return i
    def get_prolifIndex(self, compId):
        for i in range(len(self.d_proliferativeIdKeys)):
            if self.d_proliferativeIdKeys[i][:15] == compId:
                return i
    def get_mesenIndex(self, compId):
        for i in range(len(self.d_mesenchymalIdKeys)):
            if self.d_mesenchymalIdKeys[i][:15] == compId:
                return i
    def get_immunoIndex(self, compId):
        for i in range(len(self.d_immunoreactiveIdKeys)):
            if self.d_immunoreactiveIdKeys[i][:15] == compId:
                return i
                    
    def swapDashes(self, inputStr):
        dotString = ''
        for i in range(0, len(inputStr)):
            if inputStr[i] == '-':
                dotString += '.'
            else:
                dotString += inputStr[i]
        return dotString

        