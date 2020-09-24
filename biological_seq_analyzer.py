# Class for basic DNA sequence analysis

# Import files
from biological_seq_dictionaries import *

# Import modules
import random
from collections import Counter


class DNASeqAnalyzer:


    def __init__(self, sequence = 'ACGT', sequence_type = "DNA"):
        """
        Constructor for initialization of a biological sequence and the type of the sequence.
        :param sequence (str): biological sequence
        :param sequence_type (str): type of the biological sequence
        """
        self.sequence = sequence.upper()
        self.sequence_type = sequence_type

    def generate_sequence(self, nucleotides = "ACTG", sequence_length = 100):
        """
        Returns a string representing the randomly generated biological sequence.
        :param nucleotides (str):
        :param sequence_length:
        :return: str
        """
        random_sequence = "".join([random.choice(nucleotides) for i in range(sequence_length)])
        self.__init__(random_sequence)

    def sequence_length(self):
        """
        Returns an integer representing the length of a biological sequence.
        :return: int
        """
        return len(self.sequence)

    def get_seq_info(self):
        """
        Returns 2 strings representing the biological sequence and the type of the sequence.
        :return: string
        """
        return 'Sequence: {0}\nSequence type: {1}\n'.format(self.sequence, self.sequence_type)

    def nucleotide_freq_counter(self):
        """
        Returns a dictionary comprising the frequency of each nucleotide in a given sequence.
        :return: dict
        """
        for i in self.sequence:
            nucleotide_dict[i] += 1
        return nucleotide_dict

    def gc_content(self):
        """
        Returns a float representing the GC content of a biological sequence.
        :return: float
        """
        return round(((self.sequence.count('C') + self.sequence.count('G')) * 100) / len(self.sequence), 2)

    def reverse_sequence(self):
        """
        Returns a string representing the reverse version of the input biological sequence.
        :return: str
        """
        return self.sequence[::-1]

    def reverse_complement_sequence(self):
        """
        Returns a string representing the reverse complement of a DNA sequence.
        :return: str
        """
        return ''.join([dna_complement_dict[nuc] for nuc in self.sequence])[::-1]

    def dna_2_rna(self):
        """
        Returns a string representing the RNA sequence translated from a DNA sequence.
        :return: str
        """
        return self.sequence.replace('T', 'U')

    def codon_bias(self, amino_acid):
        """
        Returns a dictionary representing the frequency of codons of a given amino acid.
        :param dna_string (str): DNA string
        :param amino_acid:
        :return: dict
        """
        codon_list = []
        for i in range(0, len(self.sequence) - 2, 3):
            if dna_codon_dict[self.sequence[i:i + 3]] == amino_acid:
                codon_list.append(self.sequence[i:i + 3])

        frequency_codon_dict = dict(Counter(codon_list))
        sum_codons = sum(frequency_codon_dict.values())
        percentage_codon_dict = {codon: round(frequency_codon_dict[codon] / sum_codons, 2) for codon in
                                 frequency_codon_dict}
        return percentage_codon_dict