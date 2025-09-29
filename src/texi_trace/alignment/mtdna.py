"""
Mitochondrial DNA specific reference sequences and utilities.
"""

from typing import Dict, List, Optional

try:
    from Bio.Seq import Seq
    _BIOPYTHON_AVAILABLE = True
except ImportError:
    _BIOPYTHON_AVAILABLE = False


class MtDNAReference:
    """
    Mitochondrial DNA reference sequences and analysis utilities.
    """
    
    # Human mitochondrial DNA reference sequence (rCRS) - first 200 bases as example
    # In a full implementation, this would contain the complete 16,569 bp sequence
    HUMAN_RCRS_PARTIAL = (
        "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTG"
        "GGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTC"
        "TTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTAC"
        "TAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCAC"
    )
    
    def __init__(self, reference_sequence: Optional[str] = None):
        """
        Initialize with a reference sequence.
        
        Args:
            reference_sequence: Custom reference sequence, or None to use default
        """
        self.reference = reference_sequence or self.HUMAN_RCRS_PARTIAL
        self.length = len(self.reference)
        
    @classmethod
    def load_from_file(cls, file_path: str) -> 'MtDNAReference':
        """
        Load reference sequence from FASTA file.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            MtDNAReference instance
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        
        return cls(sequence.upper())
    
    def get_region(self, start: int, end: int) -> str:
        """
        Extract a specific region from the reference.
        
        Args:
            start: Start position (0-based)
            end: End position (0-based, exclusive)
            
        Returns:
            Sequence region
        """
        return self.reference[start:end]
    
    def get_gene_regions(self) -> Dict[str, tuple]:
        """
        Get known gene regions in human mtDNA.
        
        Returns:
            Dictionary mapping gene names to (start, end) positions
        """
        # These are approximate positions for demonstration
        # In a full implementation, these would be precise coordinates
        return {
            'D-loop': (16024, 576),     # Control region (circular)
            'rRNA_12S': (648, 1601),    # 12S ribosomal RNA
            'rRNA_16S': (1671, 3229),   # 16S ribosomal RNA  
            'ND1': (3307, 4262),        # NADH dehydrogenase subunit 1
            'COX1': (5904, 7445),       # Cytochrome c oxidase subunit I
            'COX2': (7586, 8269),       # Cytochrome c oxidase subunit II
            'ATP6': (8527, 9207),       # ATP synthase subunit 6
            'COX3': (9207, 9990),       # Cytochrome c oxidase subunit III
            'ND4': (10760, 12137),      # NADH dehydrogenase subunit 4
            'CYTB': (14747, 15887),     # Cytochrome b
        }
    
    def identify_gene_region(self, position: int) -> Optional[str]:
        """
        Identify which gene region a position falls into.
        
        Args:
            position: Position in the reference sequence
            
        Returns:
            Gene name or None if not in a known gene
        """
        gene_regions = self.get_gene_regions()
        
        for gene_name, (start, end) in gene_regions.items():
            # Handle circular genome for control region
            if gene_name == 'D-loop':
                if position >= start or position <= end:
                    return gene_name
            else:
                if start <= position <= end:
                    return gene_name
        
        return None
    
    def get_codon_usage(self, sequence: str) -> Dict[str, int]:
        """
        Analyze codon usage in a sequence.
        
        Args:
            sequence: DNA sequence to analyze
            
        Returns:
            Dictionary of codon counts
        """
        codon_counts = {}
        
        # Process sequence in codons (triplets)
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1
        
        return codon_counts
    
    def translate_sequence(self, sequence: str, genetic_code: int = 2) -> str:
        """
        Translate DNA sequence to amino acids using mitochondrial genetic code.
        
        Args:
            sequence: DNA sequence to translate
            genetic_code: NCBI genetic code table (2 = vertebrate mitochondrial)
            
        Returns:
            Amino acid sequence
        """
        # This is a simplified translation - in practice you'd use BioPython
        # with the proper mitochondrial genetic code table
        
        mitochondrial_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',  # TGA = Trp in mtDNA
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M',  # ATA = Met in mtDNA
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': '*', 'AGG': '*',  # AGA/AGG = stop in mtDNA
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        protein = ""
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                amino_acid = mitochondrial_code.get(codon, 'X')
                protein += amino_acid
        
        return protein
    
    def find_variants(self, query_sequence: str, reference_start: int = 0) -> List[Dict]:
        """
        Find variants between query and reference sequences.
        
        Args:
            query_sequence: Query sequence to compare
            reference_start: Starting position in reference
            
        Returns:
            List of variant dictionaries
        """
        variants = []
        
        end_pos = min(reference_start + len(query_sequence), len(self.reference))
        ref_segment = self.reference[reference_start:end_pos]
        
        for i, (ref_base, query_base) in enumerate(zip(ref_segment, query_sequence)):
            if ref_base != query_base:
                position = reference_start + i
                
                # Determine variant type
                variant_type = 'SNV'  # Single nucleotide variant
                
                # Check if it's in a gene region
                gene = self.identify_gene_region(position)
                
                variants.append({
                    'position': position,
                    'reference': ref_base,
                    'alternate': query_base,
                    'type': variant_type,
                    'gene': gene
                })
        
        return variants