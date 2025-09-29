class CallTrace:
    def __init__(self, trace_location, base_call, call_quality, trace_window, genome_pos=None, is_insertion=False):
        self.trace_location = trace_location  # index in trace arrays
        self.base_call = base_call            # 'A', 'C', 'G', 'T'
        self.call_quality = call_quality      # Phred score
        self.trace_window = trace_window      # list of tuples: [(A,C,G,T) at pos-5 to pos+5]
        self.genome_pos = genome_pos          # aligned genomic coordinate
        self.is_insertion = is_insertion      # flag for insertions

        self.highest_peak = max(trace_window[5]) if trace_window else 0
        self.second_highest_peak = sorted(trace_window[5], reverse=True)[1] if trace_window and len(trace_window[5]) > 1 else 0

# how can I print a CallTrace object nicely? including trace_window
    def __repr__(self):
        return (f"CallTrace(location={self.trace_location}, base={self.base_call}, "
                f"quality={self.call_quality}, genome_pos={self.genome_pos}, "
                f"is_insertion={self.is_insertion}, highest_peak={self.highest_peak}, "
                f"second_highest_peak={self.second_highest_peak}, trace_window={self.trace_window})")
