"""Sequence manipulation utilities."""

from typing import Dict


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n', '-': '-'
    }
    
    complement = ''.join(complement_map.get(base, base) for base in sequence)
    return complement[::-1]


def gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
    """
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    total_bases = len([base for base in sequence.upper() if base in 'ATGC'])
    
    if total_bases == 0:
        return 0.0
    
    return (gc_count / total_bases) * 100


def translate_dna(sequence: str, genetic_code: int = 1) -> str:
    """
    Translate DNA sequence to amino acids.
    
    Args:
        sequence: DNA sequence to translate
        genetic_code: Genetic code table (1=standard, 2=mitochondrial)
        
    Returns:
        Amino acid sequence
    """
    # Standard genetic code
    standard_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Mitochondrial genetic code differences
    mitochondrial_code = standard_code.copy()
    mitochondrial_code.update({
        'TGA': 'W',  # Tryptophan instead of stop
        'ATA': 'M',  # Methionine instead of Isoleucine
        'AGA': '*',  # Stop instead of Arginine
        'AGG': '*'   # Stop instead of Arginine
    })
    
    code_table = mitochondrial_code if genetic_code == 2 else standard_code
    
    protein = ""
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3].upper()
        if len(codon) == 3:
            amino_acid = code_table.get(codon, 'X')
            protein += amino_acid
    
    return protein


def find_orfs(sequence: str, min_length: int = 100) -> list:
    """
    Find open reading frames in a sequence.
    
    Args:
        sequence: DNA sequence to analyze
        min_length: Minimum ORF length in base pairs
        
    Returns:
        List of ORF dictionaries with start, end, and frame info
    """
    orfs = []
    
    for frame in range(3):
        for strand in [sequence, reverse_complement(sequence)]:
            strand_name = "forward" if strand == sequence else "reverse"
            
            for i in range(frame, len(strand) - 2, 3):
                codon = strand[i:i+3].upper()
                
                if codon == 'ATG':  # Start codon
                    for j in range(i + 3, len(strand) - 2, 3):
                        stop_codon = strand[j:j+3].upper()
                        
                        if stop_codon in ['TAA', 'TAG', 'TGA']:
                            orf_length = j - i + 3
                            
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': i,
                                    'end': j + 3,
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'strand': strand_name,
                                    'sequence': strand[i:j+3]
                                })
                            break
    
    return orfs