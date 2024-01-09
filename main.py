from pathlib import Path
import os

class Genbank:
    def __init__(self, FILE: Path):
        self.filename = os.path.basename(FILE)
        

class Gene:
    def __init__(self, ):
        pass


def genbank_extension_verification(FILE: Path) -> bool:
    filename = os.path.basename(FILE)
    if filename.endswith('.gbk') or filename.endswith('.gb'):
        file_is_GenBank: bool = True
        return file_is_GenBank
    else:
        raise Exception('Please provide a valid GenBank file')
    

def genbank_file_io(FILE: Path):
    file_is_GenBank = genbank_extension_verification(FILE)
    if file_is_GenBank:
        with open(FILE, 'r') as genbank_file:
            genbank_contents = genbank_file.read()       
