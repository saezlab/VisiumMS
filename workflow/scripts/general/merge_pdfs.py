import os
import fitz
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_pdfs', nargs="*", required=True)
parser.add_argument('-o','--path_merged', required=True)
args = vars(parser.parse_args())

pdfs = args['path_pdfs']
outp = args['path_merged']

result = fitz.open()
for pdf in pdfs:
    with fitz.open(pdf) as mfile:
        result.insert_pdf(mfile)
    os.remove(pdf)
result.save(outp)
