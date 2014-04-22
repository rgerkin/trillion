"""
sigma_ff_catalog.txt produced using pdf2txt.py found in PDFminer 
on sigma_ff_catalog.pdf.  
"""

with open('sigma_ff_catalog.txt','r') as f:
     text = f.read()
     lines = text.split('\n')

data = {}
organoleptic = 0
for line_num,line in enumerate(lines):
    if len(line):
        if not organoleptic and line[0]=='[':
            key = line.split(']')[0][1:]
            print key
            key = key.decode('utf-8').replace(u'\u2011','-').encode('ascii')
            organoleptic = 1
        if organoleptic and 'Organoleptic' in line:
            try:
                value = line.split(':')[1][1:]
                if value[-1] in ['-',',']:
                    if value[-1] == '-':
                        value = value[:-1]
                    else:
                        value = value + ' '
                    value += lines[line_num+1]
                value = [i.strip() for i in value.split(';') if len(i.strip())]
                data[key] = value
                organoleptic = 0
            except:
                pass

print "%d compounds described." % len(data)

descriptors = []
for x in data.values():
     descriptors += x
descriptors = list(set(descriptors)) # Remove duplicates.  
print "%d descriptors used." % len(descriptors)

