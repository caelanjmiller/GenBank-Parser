from pathlib import Path
from textwrap import wrap
import os


class GenBank:
    def __init__(self, FILE: Path):
        self.filename = os.path.basename(FILE)
        self.filepath = FILE
        self.organism: str = None
        self.accession_number: str = None
        self.genome_length: int = None
        self.genes: list = None
        self.file_contents: list = None

        def assign_organism(self, file_contents: list) -> str:
            for line in file_contents:
                if "/organism" in line:
                    self.organism = line.split("/organism=")[1].strip().replace('"', "")

        def assign_accession_number(self, file_contents: list) -> str:
            for line in file_contents:
                if "ACCESSION" in line:
                    self.accession_number = line.split("ACCESSION")[1].strip()

        def assign_genome_length(self, file_contents: list) -> int:
            for line in file_contents:
                if "REFERENCE" in line:
                    self.genome_length = line.split(" ")[-1].replace(")", "").strip()

        def annotate_cds(self, Gene, file_contents: list) -> object:
            for line in file_contents:
                if "ORIGIN" not in line:
                    raise Exception("No ORIGIN found in GenBank File")
                else:
                    pass

        with open(self.filepath, "r") as genbank_file:
            self.file_contents: list = genbank_file.readlines()
            assign_organism(self, self.file_contents)
            assign_accession_number(self, self.file_contents)
            assign_genome_length(self, self.file_contents)

    class Gene:
        def __init__(self):
            self.name: str = None
            self.coordinates = None
            self.sequence: str = None
            self.sequence_type: str = None

        def nucleotide_to_protein(self, sequence_type: str, sequence: str) -> list:
            """
            Translate nucleotide sequence to amino acid sequence
            """
            self.protein_sequence: list = None
            CODON_TABLE = {
                "GCT": "A",
                "GCC": "A",
                "GCA": "A",
                "GCG": "A",
                "TGT": "C",
                "TGC": "C",
                "GAT": "D",
                "GAC": "D",
                "GAA": "E",
                "GAG": "E",
                "TTT": "F",
                "TTC": "F",
                "GGT": "G",
                "GGC": "G",
                "GGA": "G",
                "GGG": "G",
                "CAT": "H",
                "CAC": "H",
                "ATA": "I",
                "ATT": "I",
                "ATC": "I",
                "AAA": "K",
                "AAG": "K",
                "TTA": "L",
                "TTG": "L",
                "CTT": "L",
                "CTC": "L",
                "CTA": "L",
                "CTG": "L",
                "ATG": "M",
                "AAT": "N",
                "AAC": "N",
                "CCT": "P",
                "CCC": "P",
                "CCA": "P",
                "CCG": "P",
                "CAA": "Q",
                "CAG": "Q",
                "CGT": "R",
                "CGC": "R",
                "CGA": "R",
                "CGG": "R",
                "AGA": "R",
                "AGG": "R",
                "TCT": "S",
                "TCC": "S",
                "TCA": "S",
                "TCG": "S",
                "AGT": "S",
                "AGC": "S",
                "ACT": "T",
                "ACC": "T",
                "ACA": "T",
                "ACG": "T",
                "GTT": "V",
                "GTC": "V",
                "GTA": "V",
                "GTG": "V",
                "TGG": "W",
                "TAT": "Y",
                "TAC": "Y",
                "TAA": "*",
                "TAG": "*",
                "TGA": "*",
            }
            if sequence_type == "nucleotide":
                codons = wrap(sequence, 3)
                for codon in codons:
                    if codon in CODON_TABLE.keys():
                        self.protein_sequence.append(CODON_TABLE[codon])
