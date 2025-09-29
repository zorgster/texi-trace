class AnalysisError(Exception):
    """Raised when an analysis step fails."""
    def __init__(self, message, file=None):
        super().__init__(message)
        self.file = file

# raise AnalysisError("FASTQ conversion failed", file=abi_file)

# except AnalysisError as e:
#    print(f"[ERROR] {e.file}: {e}", file=sys.stderr)
