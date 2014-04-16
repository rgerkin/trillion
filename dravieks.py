"""
dravieks.txt produced by copying and pasting from the supplemental material of 
"In search of the structure of human olfactory space", e.g. koulakov.pdf.  
"""

with open('dravieks.txt','r') as f:
    text = f.read()
    lines = text.split('\n')
    dravieks = [line.split(' ')[1] for line in lines] # A list of CAS numbers.  
    