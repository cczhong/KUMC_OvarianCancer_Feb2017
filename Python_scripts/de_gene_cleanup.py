import sys


    
def geneNameCleanup(file):
    geneName_output = []
    with open(file,  'r') as f:
        read_str = ''
        
        for line in f:
            read_str = line.rstrip()
            i = 0
            clean_str= ''
            pipe_found = False
            
            while pipe_found == False:
                if read_str[i] == '|':
                    geneName_output.append(clean_str)
                    pipe_found = True
                    
                else:
                    if read_str[i] != '"':
                        clean_str += read_str[i]
                i += 1
              
    return geneName_output

f_names = geneNameCleanup(sys.argv[1])
m_names = geneNameCleanup(sys.argv[2])
i_names = geneNameCleanup(sys.argv[3])
p_names = geneNameCleanup(sys.argv[4])

f_out = open('fallopianDE_names.txt', 'w+')
m_out = open('mesenchymalDE_names.txt', 'w+')
i_out = open('immunoreactiveDE_names.txt','w+')
p_out = open('proliferativeDE_names.txt', 'w+')

for i in range(0, len(f_names)):
    f_out.write(f_names[i])
    f_out.write('\n')
for i in range(0, len(m_names)):
    m_out.write(m_names[i])
    m_out.write('\n')
for i in range(0, len(i_names)):
    i_out.write(i_names[i])
    i_out.write('\n')
for i in range(0, len(p_names)):
    p_out.write(p_names[i])
    p_out.write('\n')    

f_out.close()
m_out.close()
i_out.close()
p_out.close()
