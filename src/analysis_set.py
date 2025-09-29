import pickle
from datetime import datetime
from pathlib import Path
import getpass

class AnalysisSet:
    def __init__(self, name, main_ratio=0.75, secondary_ratio=0.3, quality_threshold=20, start =None, end=None):
        self.entries = []
        self.main_ratio = main_ratio
        self.secondary_ratio = secondary_ratio
        self.quality_threshold = quality_threshold
        self.start = start
        self.end = end

        self.metadata = {
            'created_at': datetime.now(),
            'user': getpass.getuser(),
            'version': 1.0,
        }

        self.name = name

        #Counters used along the script
        self.analysed = 0
        self.possibilities = 0
        self.lowquality = 0
        self.Averagequality = list()


    def add_file(self, file_path):
        self.files.append(file_path)

    def remove_file(self, file_path):
        self.files.remove(file_path)

    def list_files(self):
        return self.files

    def __repr__(self):
        return f"AnalysisSet(name={self.name}, files={self.files})"

