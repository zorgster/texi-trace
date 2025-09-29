from analysis_error import AnalysisError
from collections import namedtuple
import re
import subprocess
from os import path
from Bio import SeqIO

SoftClipOffset = namedtuple('SoftClipOffset', ['start', 'end'])

class SamAlignment:
    def __init__(self, sam_file, index_path, force_overwrite=False):
        self.sam_file = sam_file
        self.fastq_file = sam_file.replace(".sam", ".fastq")

        self.index_path = index_path
        self.force_overwrite = force_overwrite

        # Alignment details
        self.ref_name = None  # Reference sequence name
        self.cigar = [] # List of (length, operation) tuples
        self.is_forward = True # True if forward strand
        self.reference_position = 1  # Position in reference (1-based)
        self.mapping_quality = 0
        self.softclip_offset = SoftClipOffset(start=0, end=0) # Number of clipped bases from SAM CIGAR

        self.genome_map = {}  # Maps read index (0-based) to genome position (1-based)
        self.reverse_genome_map = {} # Maps genome position (1-based) to read index (0-based)

    def _process(self):
        self._align_fastq_to_index()
        self._parse_sam()
        self._map_read_to_genome()

        self.valid = True if self.cigar else False
        return self

    @classmethod
    def create(cls, sam_file, index_path, force_overwrite=False):
        sam = cls(sam_file=sam_file, index_path=index_path, force_overwrite=force_overwrite)

        try:
            sam._process()
            if sam.valid:
                return sam
            else:
                return AnalysisError("Invalid SAM alignment.", file=sam_file)
        except AnalysisError as e:
            return AnalysisError(f"SAM processing failed: {e}", file=sam_file)

    # Align FASTQ to reference index using Bowtie2 to produce SAM file
    def _align_fastq_to_index(self):
        # if SAM already exists and not forcing, skip
        if path.exists(self.sam_file) and not self.force_overwrite:
            return

        # check index exists
        required_index_files = [f"{self.index_path}.{i}.bt2" for i in range(1, 5)]
        for f in required_index_files:
            if not path.exists(f):
                raise AnalysisError(f"Missing reference index: {f}", file=self.sam_file)

        # Call method to convert AB1 to FASTQ (if required)
        self._convert_ab1_to_fastq()

        # Call Bowtie2 to align FASTQ to index, output SAM
        bowtie2_cmd = ["bowtie2", "--local",
                        "--mp", "9,5",  # Mismatch penalty
                        "--rdg", "7,3", # Read gap open, extend 
                        "--rfg", "7,3", # Ref gap open, extend
                        "-q", self.fastq_file, 
                        "-x", self.index_path, 
                        "-S", self.sam_file]
        try:
            subprocess.run(bowtie2_cmd, check=True)
        except Exception as e:
            raise AnalysisError(f"Bowtie2 alignment failed: {e}", file=self.sam_file)

    # Convert AB1 to FASTQ using Biopython
    def _convert_ab1_to_fastq(self):
        # Check if FASTQ already exists and not forcing
        if path.exists(self.fastq_file) and not self.force_overwrite:
            return True

        ab1_file = self.sam_file.replace(".sam", ".ab1")

        try:
            SeqIO.convert(ab1_file, "abi", self.fastq_file, "fastq")
            return True
        except Exception as e:
            raise AnalysisError(f"Failed to convert AB1 to FASTQ: {e}", file=self.sam_file)

    # Parse SAM file to extract alignment details
    def _parse_sam(self):
        """Parse SAM file from single AB1 read."""
        # Check again that SAM file exists?  It should have been created earlier. Call _align_fastq_to_index() first
        if not path.exists(self.sam_file):
            raise AnalysisError(f"SAM file missing: {self.sam_file}", file=self.sam_file)

        # Read SAM file, extract alignment info from first non-header line
        try:
            line_processed = False
            with open(self.sam_file, 'r') as f:
                for line in f:
                    if not line.startswith('@'):
                        fields = line.strip().split('\t')
                        if len(fields) > 5:
                            # keep strand as int
                            self.is_forward = (fields[1] == '0') # True if forward strand
                            self.ref_name = fields[2]  # Reference sequence name
                            self.reference_position = int(fields[3])  # Position in reference (1-based)
                            self.mapping_quality = int(fields[4])
                            # keep softclip_offset as int (0-based) ('extra' in old version)
                            self.cigar = self._parse_cigar(fields[5])
                            self.softclip_offset = self._parse_softclip_offset(self.cigar)
                            line_processed = True
                        else:
                            raise AnalysisError(f"Malformed SAM entry: {line.strip()}", file=self.sam_file)
            if not line_processed:
                raise AnalysisError(f"No alignment entries found in SAM file", file=self.sam_file)
        except Exception as e:
            raise AnalysisError(f"Failed to parse SAM file - {e}", file=self.sam_file)

    def _parse_cigar(self, cigar_str: str) -> list:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)]
    
    def _parse_softclip_offset(self, cigar: list) -> int:
        """Extract leading and trailing softclip length from CIGAR list."""
        start_match = cigar[0][0] if cigar and cigar[0][1] == 'S' else 0
        end_match = cigar[-1][0] if cigar and cigar[-1][1] == 'S' else 0
        return SoftClipOffset(start=start_match, end=end_match)

    def _map_read_to_genome(self):
        read_to_genome = {}
        read_pos = 0
        genome_pos = self.reference_position

        # The CIGAR and reference_position are always given in the 5' to 3' direction of the reference
        # The read positions are always 0-based from the start of the read, regardless of strand
        # So we can just iterate through the CIGAR and map read positions to genome positions directly

        # It is the responsibility of the caller to handle strand direction if needed

        for length, op in self.cigar:
            if op in ('M', '=', 'X'):
                for _ in range(length):
                    read_to_genome[read_pos] = genome_pos
                    read_pos += 1
                    genome_pos += 1
            elif op == 'I':
                for _ in range(length):
                    read_to_genome[read_pos] = genome_pos  # insertion: same genome_pos
                    read_pos += 1
            elif op == 'D':
                genome_pos += length  # deletion: skip reference bases
            elif op == 'S':
                read_pos += length  # soft clip: skip read bases
            elif op == 'H':
                continue  # hard clip: not in read

        self.genome_map = read_to_genome

        # Create reverse map too
        self.reverse_genome_map = {v: k for k, v in read_to_genome.items()}


    def get_genome_position(self, read_index):
        """Get genomic position for a given read index (0-based)."""
        return self.genome_map.get(read_index, None)

    # Find read_index for a genome position
    def get_read_index(self, genome_position):
        """Get read index (0-based) for a given genomic position (1-based)."""
        return self.reverse_genome_map.get(genome_position, None)
