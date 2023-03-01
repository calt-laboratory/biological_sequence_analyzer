##### Class for basic DNA sequence analysis
##### Author: Christoph Alt

# Import modules
import random
from collections import Counter
from typing import Dict, Union

import matplotlib.pyplot as plt

# Import files
from biological_seq_dictionaries import *

class BioSeqAnalyzer:
    """
    Class for analyzing biological sequences (RNA, DNA).
    """

    def __init__(self, sequence: str = "ACGT", sequence_type: str = "DNA") -> None:
        """
        Constructor for initialization of a biological sequence and the type of the sequence.
        """
        self.sequence = sequence.upper()
        self.sequence_type = sequence_type

    def __repr__(self) -> None:
        """
        Returns a formal representation of an instanced object as string.
        """
        print(f'BiolgicalSeqAnalyzer(sequence={self.sequence}, sequence_type={self.sequence_type}')

    def __str__(self) -> None:
        """
        Returns an informal representation of an instanced as string.
        """
        print(f'Biolgical Sequence Analyzer: Seq = {self.sequence}, Seq_Type = {self.sequence_type}')

    def generate_sequence(self, nucleotides: str = "ACGT", sequence_type: str = "DNA", sequence_length: str = 100):
        """
        Returns a string representing the randomly generated biological sequence.
        """
        random_sequence = "".join([random.choice(nucleotides) for i in range(sequence_length)])
        self.__init__(random_sequence, sequence_type)

    def sequence_length(self) -> int:
        """
        Returns an integer representing the length of a biological sequence.
        """
        return len(self.sequence)

    def get_seq_info(self) -> str:
        """
        Returns 2 strings representing the biological sequence and the type of the sequence.
        """
        return 'Sequence: {0}\nSequence type: {1}\n'.format(self.sequence, self.sequence_type)

    def nucleotide_freq_counter(self) -> Dict[str, int]:
        """
        Returns a dictionary comprising the frequency of each nucleotide in a given sequence.
        """
        if self.sequence_type == "DNA":
            for i in self.sequence:
                dna_nucleotide_dict[i] += 1
            return dna_nucleotide_dict
        if self.sequence_type == "RNA":
            for i in self.sequence:
                rna_nucleotide_dict[i] += 1
            return rna_nucleotide_dict

    def gc_content(self) -> float:
        """
        Returns a float representing the GC content of a biological sequence.
        """
        return round(((self.sequence.count('C') + self.sequence.count('G')) * 100) / len(self.sequence), 2)

    def reverse_sequence(self) -> str:
        """
        Returns a string representing the reverse version of the input biological sequence.
        """
        return self.sequence[::-1]

    def reverse_complement_sequence(self) -> str:
        """
        Returns a string representing the reverse complement of a biological sequence.
        """
        if self.sequence_type == "DNA":
            return ''.join([dna_complement_dict[nuc] for nuc in self.sequence])[::-1]
        elif self.sequence_type == "RNA":
            return ''.join([rna_complement_dict[nuc] for nuc in self.sequence])[::-1]

    def dna_2_rna(self) -> str:
        """
        Returns a string representing the RNA sequence translated from a DNA sequence.
        """
        if self.sequence_type == "DNA":
            return self.sequence.replace('T', 'U')
        else:
            return "Given sequence is not a DNA sequence!"

    def codon_bias(self, amino_acid: str = "A", viz: bool = True) -> Dict:
        """
        Returns a dictionary representing the frequency of codons of a given amino acid.
        """
        codon_list = []
        if self.sequence_type == "DNA":
            for i in range(0, len(self.sequence) - 2, 3):
                if dna_codon_dict[self.sequence[i:i + 3]] == amino_acid:
                    codon_list.append(self.sequence[i:i + 3])
        if self.sequence_type == "RNA":
            for i in range(0, len(self.sequence) - 2, 3):
                if rna_codon_dict[self.sequence[i:i + 3]] == amino_acid:
                    codon_list.append(self.sequence[i:i + 3])

        frequency_codon_dict = dict(Counter(codon_list))
        sum_codons = sum(frequency_codon_dict.values())
        percentage_codon_dict = {codon: round(frequency_codon_dict[codon] / sum_codons, 2) for codon in
                                 frequency_codon_dict}

        if viz == True:
            labels = list(percentage_codon_dict.keys())
            vals = list(percentage_codon_dict.values())
            plt.pie(vals, labels=vals , normalize=False)
            plt.legend(labels)
            plt.title(f'Frequency of codons for the amino acid: {amino_acid}')
            plt.show()
        return percentage_codon_dict

    def rna_2_protein(self) -> Union[str, None]:
        """
        Returns a string representing a amino acid sequence translated from a RNA sequence.
        """
        if self.sequence_type == "RNA":
            print(type(''.join([rna_codon_dict[self.sequence[i:i + 3]] for i in range(0, len(self.sequence) - 2, 3)])))
            return ''.join([rna_codon_dict[self.sequence[i:i + 3]] for i in range(0, len(self.sequence) - 2, 3)])
        else:
            print('Given sequence is not a RNA!')

def main():
    seq_analyzer = BioSeqAnalyzer(sequence="ACUGAUUUUUACCCAAAA", sequence_type="RNA")
    seq_analyzer.rna_2_protein()


if __name__ == "__main__":
    main()
