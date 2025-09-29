"""File handling utilities."""

import os
from pathlib import Path
from typing import List, Optional


def validate_file_format(file_path: str, allowed_extensions: List[str]) -> bool:
    """
    Validate that a file has an allowed extension.
    
    Args:
        file_path: Path to the file
        allowed_extensions: List of allowed extensions (e.g., ['.ab1', '.scf'])
        
    Returns:
        True if file extension is allowed
    """
    ext = Path(file_path).suffix.lower()
    return ext in [e.lower() for e in allowed_extensions]


def create_output_directory(output_path: str) -> str:
    """
    Create output directory if it doesn't exist.
    
    Args:
        output_path: Path to output directory
        
    Returns:
        Absolute path to created directory
    """
    Path(output_path).mkdir(parents=True, exist_ok=True)
    return os.path.abspath(output_path)


def find_files_by_extension(directory: str, extensions: List[str]) -> List[str]:
    """
    Find all files with specified extensions in a directory.
    
    Args:
        directory: Directory to search
        extensions: List of file extensions to find
        
    Returns:
        List of file paths
    """
    directory_path = Path(directory)
    files = []
    
    for ext in extensions:
        pattern = f"*{ext}"
        files.extend(directory_path.glob(pattern))
    
    return [str(f) for f in files]


def get_file_size(file_path: str) -> int:
    """
    Get file size in bytes.
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in bytes
    """
    return os.path.getsize(file_path)