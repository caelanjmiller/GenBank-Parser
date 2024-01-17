from pathlib import Path
from textwrap import wrap
from itertools import islice
import os


class GenBank:
    def __init__(self, FILE: Path):
        self.filename = os.path.basename(FILE)
        self.filepath = FILE
        self.organism: str = None
        self.accession_number: str = None
        self.genome_length: int = None
        self.genes: list = []
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
                    self.genome_length = int(line.split(" ")[-1].replace(")", "").strip())

        def annotate_cds(self, file_contents: list) -> object:
            # Capture indices at which gene feature key occurs in GenBank file content list
            gene_feature_indices: list = [index for index, value in enumerate(file_contents) if value.startswith('gene')]
            # Capture index at which ORIGIN occurs - signals end of GenBank annotation before full sequence
            origin_index: int = [index for index, value in enumerate(file_contents) if value.startswith('ORIGIN')][0]
            # Create a list containing tuples of the indices of the information contained between each occurrence of gene feature keys
            information_indices: list = [tuple((gene_feature_indices[index] + 1, gene_feature_indices[index + 1])) for index in range(len(gene_feature_indices) - 1)]
            # Capture last gene entry prior to ORIGIN
            information_indices.append(tuple((gene_feature_indices[-1], origin_index)))
            parsed_file_contents: list = [list(islice(file_contents, a, b)) for a, b in information_indices]
            for gene_slice in parsed_file_contents:
                gene = Gene()
                gene.assign_annotated_function(gene_slice)
                gene.assign_protein_id(gene_slice)
                gene.assign_sequence(gene_slice)
                self.genes.append(gene)


        with open(self.filepath, "r") as genbank_file:
            self.file_contents: list = [line.strip() for line in genbank_file.readlines()]
            assign_organism(self, self.file_contents)
            assign_accession_number(self, self.file_contents)
            assign_genome_length(self, self.file_contents)
            annotate_cds(self, self.file_contents)

class Gene:
    def __init__(self):
        self.protein_id: str = None
        self.annotated_function: str = None
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
    
    def assign_annotated_function(self, file_contents: list) -> str:
        for line in file_contents:
            if '/product' in line:
                self.annotated_function = line.split('/product=')[1].strip().replace('"', "")
    
    def assign_protein_id(self, file_contents: list) -> str:
        for line in file_contents:
            if '/protein_id' in line:
                self.protein_id = line.split('/protein_id=')[1].strip().replace('"', "")
            else:
                self.protein_id = 'NA'

    def assign_coordinates(self, file_contents: list) -> tuple:
        pass

    def assign_sequence(self, file_contents: list) -> str:
        for gene_entry in file_contents:
            if '/translation=' in gene_entry:
                translation_feature_index: int = [index for index, value in enumerate(file_contents) if ('/translation=') in value][0]
                end_sequence_indices: list = ([index for index, value in enumerate(file_contents) if re.search(r'[A-Z]"\Z', value)])
                lines_after_translation: list = []
                if len(end_sequence_indices) > 1:
                    end_sequence_index: int = end_sequence_indices[1] + 1
                else:
                    end_sequence_index: int = end_sequence_indices[0] + 1
            for line in file_contents[translation_feature_index:end_sequence_index]:
                lines_after_translation.append(line)
                # Join list of translation sequence together into string and format by removing whitespace and quotations
                sequence: str = "".join(lines_after_translation).split('/translation=')[1].strip().replace('"', "")
                self.sequence = sequence
            else:
                continue