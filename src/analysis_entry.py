from datetime import datetime
from os import path
from analysis_error import AnalysisError
from call_trace import CallTrace
from Bio import SeqIO
from Bio.Data import IUPACData
import re
from sys import stderr
import pickle
from collections import namedtuple
import numpy as np
import math
from sam_alignment import SamAlignment
import matplotlib.pyplot as plt

class AnalysisEntry:
    ROI = namedtuple('ROI', ['start', 'end'])
    PeakNeighbours = namedtuple('PeakNeighbours', ['before1', 'before2', 'before3', 'after1', 'after2', 'after3'])

    def __init__(self, ab1_file, start=100, end=150, quality_threshold=20):
        # File paths
        self.ab1_file = ab1_file
        self.sam_file = self.ab1_file.replace(".ab1", ".sam")

        # Analysis parameters
        self.roi = self.ROI(start=start, end=end)
        self.mean_quality_threshold = quality_threshold

        self.call_traces = []  # List of CallTrace objects

        # Trace / basecall data
        self.trace_A = []
        self.trace_T = []
        self.trace_G = []
        self.trace_C = []
        self.quality_scores = []
        self.base_calls = []
        self.call_locations = []
        self.count_of_calls = 0 # Length of the sequence in the AB1 file
        self.trace_at_calls = []  # List of PeakNeighbours namedtuples at call positions
        self.mean_quality = 0 # Average quality score across all base calls in aligned read
        self.mean_quality_in_roi = 0 # Average quality score across base calls in region of interest (start to end)

        # Alignment info
        self.sam_alignment = None # SamAlignment object

        # Status flags
        self.valid = False
        self.error_log = []

    def _process(self, index_path, force_overwrite=False):
        self.sam_alignment = SamAlignment.create(self.sam_file, index_path, force_overwrite)
        if isinstance(self.sam_alignment, AnalysisError):
            raise self.sam_alignment

        self._extract_trace_data()
        self._extract_trace_at_calls()
        self._extract_call_traces()
        self._extract_quality_at_calls()
        # TODO: Filtering possibilities based on quality and alignment
        self._filter_possibilities()
        stderr.write(f"Mean quality in ROI (start={self.roi.start}, end={self.roi.end}): {self.mean_quality_in_roi}\n")
        self.valid = True
        return self

    @classmethod
    def create(cls, ab1_file, ref_index, force=False, start=100, end=150, quality_threshold=20):
        # Verify AB1 file exists
        if not path.exists(ab1_file):
            print(f"[ERROR] {ab1_file}: AB1 file missing", file=stderr)
            return None

        entry = cls(ab1_file, start=start, end=end, quality_threshold=quality_threshold)

        try:
            entry._process(index_path=ref_index, force_overwrite=force)
            stderr.write(f"Processed {ab1_file} successfully.\n")
            return entry if entry.valid else None
        except AnalysisError as e:
            print(f"[ERROR] {e.file}: {e}", file=stderr)
            return None

    # Pickle methods
    def save(self, pickle_file: str = None):
        ''' Save AnalysisEntry to pickle file '''
        if not pickle_file:
            pickle_file = self.ab1_file.replace(".ab1", ".entry")
        with open(pickle_file, 'wb') as f:
            pickle.dump(self, f)
        print(f"Saved AnalysisEntry to {pickle_file}")

    @staticmethod
    def load(pickle_file: str):
        ''' Load AnalysisEntry from pickle file '''
        with open(pickle_file, 'rb') as f:
            return pickle.load(f)

    def _extract_trace_data(self):
        """Extract trace and basecall data from AB1 file."""
        try:
            record = SeqIO.read(self.ab1_file, "abi")
            abif = record.annotations.get('abif_raw', {})

            if not abif:
                raise AnalysisError(f"No ABIF data found in {self.ab1_file}", file=self.ab1_file)

            # Extract trace data
            self.trace_G = abif.get('DATA9', [])
            self.trace_A = abif.get('DATA10', [])
            self.trace_T = abif.get('DATA11', [])
            self.trace_C = abif.get('DATA12', [])

            self.quality_scores = list(abif.get('PCON2', []))
            self.base_calls = abif.get('PBAS2', [])
            self.call_locations = abif.get('PLOC2', [])
            self.count_of_calls = len(self.base_calls)

            del record  # Free memory
            del abif  # Free memory

            if not (self.trace_A and self.trace_T and self.trace_G and self.trace_C):
                raise AnalysisError(f"Incomplete trace data in {self.ab1_file}", file=self.ab1_file)

            if not (len(self.trace_A) == len(self.trace_T) == len(self.trace_G) == len(self.trace_C)):
                raise AnalysisError(f"Mismatched trace data lengths in {self.ab1_file}", file=self.ab1_file)

            if not (self.quality_scores and self.base_calls and self.call_locations):
                raise AnalysisError(f"Incomplete basecall data in {self.ab1_file}", file=self.ab1_file)

            if not (len(self.quality_scores) == len(self.base_calls) == len(self.call_locations)):
                raise AnalysisError(f"Mismatched basecall data lengths in {self.ab1_file}", file=self.ab1_file)

            return True
        except Exception as e:
            raise AnalysisError(f"Failed to extract trace data: {e}", file=self.ab1_file)

    def _extract_trace_at_calls(self):
        """Extract trace values at base call positions."""
        if not (self.trace_A and self.trace_T and self.trace_G and self.trace_C):
            raise AnalysisError(f"Trace data not available for {self.ab1_file}", file=self.ab1_file)

        # Get clipped calls - returns the forward or reverse set of call locations removing the initial softclip
        clipped_calls = self._get_clipped_calls()

        # Get tuples of trace values at each call position
        self.trace_at_calls = [
            (self.trace_A[i], self.trace_C[i], self.trace_G[i], self.trace_T[i]) for i in clipped_calls
        ]

    def _extract_call_traces(self):
        """Create CallTrace objects for each base call."""
        if not self.sam_alignment:
            raise AnalysisError(f"SAM alignment not available for {self.ab1_file}", file=self.ab1_file)

        # Get clipped calls (locations) - returns the forward or reverse set of call locations removing the softclipped bases
        clipped_calls = self._get_clipped_calls()
        # Get clipped bases - returns the forward or reverse complement of base calls removing the softclipped bases
        clipped_bases = self._get_clipped_bases()
        # Get clipped qualities - returns the forward or reverse set of quality scores removing the softclipped bases
        clipped_qualities = self._get_clipped_qualities()
        
        # Fill call_traces list for every called base
        for i in range(len(clipped_calls)):
                call_trace = CallTrace(
                    trace_location=clipped_calls[i],
                    base_call=clipped_bases[i],
                    call_quality=clipped_qualities[i],
                    trace_window=self._safe_trace_window(clipped_calls[i], flank=5),
                    genome_pos=self.sam_alignment.get_genome_position(i) if self.sam_alignment else -1,
                    is_insertion=False  # Placeholder; actual insertion detection not implemented yet
                )
                self.call_traces.append(call_trace)

    def _get_clipped_calls(self):
        """Return list of call locations adjusted for strand and softclip."""
        if self.sam_alignment.is_forward:
            return self.call_locations[self.sam_alignment.softclip_offset.start:self.count_of_calls - self.sam_alignment.softclip_offset.end]
        else:
            # For reverse strand, reverse the call locations and adjust
            clipped = self.call_locations[self.sam_alignment.softclip_offset.end:self.count_of_calls - self.sam_alignment.softclip_offset.start]
            return list(reversed(clipped))

    def _get_clipped_bases(self):
        """Return list of base calls adjusted for strand and softclip."""
        if self.sam_alignment.is_forward:
            return [chr(b) for b in self.base_calls[self.sam_alignment.softclip_offset.start:self.count_of_calls - self.sam_alignment.softclip_offset.end]]
        else:
            # For reverse strand, reverse the base calls and complement
            clipped = self.base_calls[self.sam_alignment.softclip_offset.end:self.count_of_calls - self.sam_alignment.softclip_offset.start]
            IUPAC = IUPACData.ambiguous_dna_complement
            return [IUPAC[chr(b)] for b in reversed(clipped)]

    def _get_clipped_qualities(self):
        """Return list of quality scores adjusted for strand and softclip."""
        if self.sam_alignment.is_forward:
            return self.quality_scores[self.sam_alignment.softclip_offset.start:self.count_of_calls - self.sam_alignment.softclip_offset.end]
        else:
            # For reverse strand, reverse the quality scores
            clipped = self.quality_scores[self.sam_alignment.softclip_offset.end:self.count_of_calls - self.sam_alignment.softclip_offset.start]
            return list(reversed(clipped))

    def _extract_quality_at_calls(self):
        """Extract quality scores at base call positions."""
        if not self.quality_scores:
            raise AnalysisError(f"Quality scores not available for {self.ab1_file}", file=self.ab1_file)

        clipped_quals = self._get_clipped_qualities()

        if clipped_quals:
            mean_p = np.mean([10 ** (-q / 10) for q in clipped_quals])
            self.mean_quality = -10 * math.log10(mean_p)

            # Adjust for 0-based indexing
            start_idx = self.roi.start - 1
            end_idx = self.roi.end

            # Filter qualities to those within the ROI
            roi_quals = [q for i, q in enumerate(clipped_quals) if start_idx <= i < end_idx]
            stderr.write(f"ROI Quals: {roi_quals}\n")

            # Calculate mean quality in ROI
            if roi_quals:
                mean_p = np.mean([10 ** (-q / 10) for q in roi_quals])
                self.mean_quality_in_roi = -10 * math.log10(mean_p)
            else:
                self.mean_quality_in_roi = 0

        else:
            self.mean_quality = 0
            self.mean_quality_in_roi = 0

        # No need to add this to Averagequality list in AnalysisSet; that should be handled externally

    def _filter_possibilities(self):
        """Filter possibilities based on quality and alignment criteria."""
        # Placeholder for actual filtering logic
        # For now, just set possibilities to 1 if mean quality in ROI exceeds threshold
        if self.mean_quality_in_roi < self.mean_quality_threshold:  # Example threshold
            raise AnalysisError(f"Low quality in ROI: {self.mean_quality_in_roi}", file=self.ab1_file)
        
        # Additional filtering based on alignment can be added here

        for i, (peak) in enumerate(self.trace_at_calls):
            # Example filtering logic: count peaks above a certain threshold
            if i in range(self.roi.start, self.roi.end):
                # TODO: try:
                    #     peaks = PeakNeighbours(before1=peak[i-1] if i-1 >= 0 else 0,
                    #                            before2=peak[i-2] if i-2 >= 0 else 0,
                    #                            before3=peak[i-3] if i-3 >= 0 else 0,
                    #                            after1=peak[i+1] if i+1 < len(peak) else 0,
                    #                            after2=peak[i+2] if i+2 < len(peak) else 0,
                    #                            after3=peak[i+3] if i+3 < len(peak) else 0)
                # except Exception as e:
                #     raise AnalysisError(f"Error extracting peak neighbours: {e}", file=self.ab1_file)

                # TODO: Get prior and post peaks
                # TODO: Sort CurrentPeak - 1. highest is [3], 2. lower = [2]
                # TODO: Ratio = (lower / highest) * 100
                # TODO: get Lowerb1 = PeakHeight(before1) LowerPeak
                # TODO: get Lowerb2 = PeakHeight(before2) LowerPeak
                # TODO: get Lowerb3 = PeakHeight(before3) LowerPeak
                # TODO: get Lowera1 = PeakHeight(after1) LowerPeak
                # TODO: get Lowera2 = PeakHeight(after2) LowerPeak
                # TODO: get Lowera3 = PeakHeight(after3) LowerPeak

                # TODO: Test conditions for heteroplasmy to consider double peak
                # TODO: Ratio >= MainRatio AND Lowerb1 <= (Lower * self.Peakb1)
                #      AND Lowerb2 <= (Lower * self.Peakb2)
                #      AND Lowerb3 <= (Lower * self.Peakb3)
                #      AND Lowera1 <= (Lower * self.Peaka1)
                #      AND Lowera2 <= (Lower * self.Peaka2)
                #      AND Lowera3 <= (Lower * self.Peaka3):
                #          position = i + 1 + inicio? ; counter += 1 possibilities += 1
                # Funciona?
                # if self.is_forward:
                #    hp.write - chroma, position, Ratio, MeanQuality\n
                # or reverse:
                #    hp.write - ... same tab Reverse\n

                return

    def _safe_trace_window(self, i, flank=5):
        start = max(0, i - flank)
        end = min(len(self.trace_A), i + flank + 1)
        return [
            (self.trace_A[j], self.trace_C[j], self.trace_G[j], self.trace_T[j])
            for j in range(start, end)
        ]

    def is_heteroplasmic(self, index, threshold=0.3):
        signals = {
            'A': self.trace_A[index],
            'C': self.trace_C[index],
            'G': self.trace_G[index],
            'T': self.trace_T[index]
        }
        sorted_peaks = sorted(signals.items(), key=lambda x: x[1], reverse=True)
        stderr.write(f"Index {index} peaks: {sorted_peaks}\n")

        top, second = sorted_peaks[0], sorted_peaks[1]
        stderr.write(f"Top peak: {top}, Second peak: {second}\n")

        # Ratio of second peak to top peak
        ratio = second[1] / top[1] if top[1] > 0 else 0

        return ratio >= threshold, top[0], second[0], ratio

    def plot_at_positions(self, call_location, width=5, filename="trace_plot.png"):
        # Define the trace indices and base calls
        center_indices = [self.call_locations[i] for i in range(call_location-width, call_location+width+1) if 0 <= i < len(self.call_locations)]
        base_calls = [chr(self.base_calls[i]) for i in range(call_location-width, call_location+width+1) if 0 <= i < len(self.base_calls)]

        # Expand each index ±5
        expanded_indices = sorted(set(i for idx in center_indices for i in range(idx - 5, idx + 6)))

        # Extract trace values
        trace_A = [self.trace_A[i] for i in expanded_indices]
        trace_C = [self.trace_C[i] for i in expanded_indices]
        trace_G = [self.trace_G[i] for i in expanded_indices]
        trace_T = [self.trace_T[i] for i in expanded_indices]

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(expanded_indices, trace_A, label='A', color='green', marker='o')
        plt.plot(expanded_indices, trace_C, label='C', color='blue', marker='o')
        plt.plot(expanded_indices, trace_G, label='G', color='black', marker='o')
        plt.plot(expanded_indices, trace_T, label='T', color='red', marker='o')

        # Annotate base calls at center positions
        for i, idx, base in zip(
            range(call_location - width, call_location + width + 1),
            center_indices,
            base_calls
        ):
            if 0 <= i < len(self.call_locations):
                gpos = self.sam_alignment.get_genome_position(i)
                y = max(self.trace_A[idx], self.trace_C[idx], self.trace_G[idx], self.trace_T[idx])
                label = f"{base}\n{gpos}"
                plt.text(idx, y + 50, label, ha='center', va='bottom', fontsize=9, fontweight='bold')

        plt.title("Trace Signal Across Genomic Region")
        plt.xlabel("Trace Index")
        plt.ylabel("Signal Intensity")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        #plt.show()

