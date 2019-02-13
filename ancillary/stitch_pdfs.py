import time
import sys
import os
import subprocess as sp
import numpy as np

write_path = sys.argv[1]

os.chdir(write_path)

#get list of PDFs to stitch
dirfiles = os.listdir(write_path)
indiv_pdfs = list(filter(lambda x: '.pdf' in x, dirfiles)) #filter non-pdfs

#sort by filename
indiv_pdfs = list(np.array(indiv_pdfs)[np.argsort(indiv_pdfs)])

#execute `echo <files> out.pdf | xargs pdfunite`
p1 = sp.Popen(['echo'] + indiv_pdfs + ['out.pdf'], stdout=sp.PIPE)
p1.wait() #wait for subprocess to finish before executing next one
p2 = sp.Popen(['xargs', 'pdfunite'], stdin=p1.stdout)
p2.wait()

