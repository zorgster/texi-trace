"""
Command-line interface for texi-trace.
"""

import click
import os
from pathlib import Path
from typing import Optional

from .chromatogram import ChromatogramReader
from .alignment import SequenceAligner
from .visualization import ChromatogramVisualizer


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """
    texi-trace: A tool for scrutinizing Sanger chromatograms for mtDNA analysis.
    
    Provides functionality to analyze chromatogram files, perform sequence alignment
    to reference genomes, and visualize results.
    """
    pass


@cli.command()
@click.argument('chromatogram_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', type=click.Path(), default='output',
              help='Output directory for results')
@click.option('--format', '-f', type=click.Choice(['abi', 'scf']), 
              help='Chromatogram file format (auto-detected if not specified)')
def analyze(chromatogram_file: str, output_dir: str, format: Optional[str]):
    """
    Analyze a chromatogram file and generate basic reports.
    
    CHROMATOGRAM_FILE: Path to the chromatogram file (.ab1, .abi, .scf)
    """
    click.echo(f"Analyzing chromatogram file: {chromatogram_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read chromatogram
        reader = ChromatogramReader(chromatogram_file)
        data = reader.read()
        
        # Basic info
        click.echo(f"Sequence length: {len(data['sequence'])} bases")
        click.echo(f"Average quality: {sum(data['quality_scores'])/len(data['quality_scores']):.2f}")
        
        # Generate visualization
        visualizer = ChromatogramVisualizer()
        
        # Plot traces
        if data['traces']:
            fig = visualizer.plot_traces(data['traces'], data['sequence'])
            output_path = os.path.join(output_dir, 'traces.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            click.echo(f"Trace plot saved to: {output_path}")
        
        # Plot quality scores
        if data['quality_scores']:
            fig = visualizer.plot_quality_scores(data['quality_scores'], data['sequence'])
            output_path = os.path.join(output_dir, 'quality.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            click.echo(f"Quality plot saved to: {output_path}")
        
        # Save sequence to file
        seq_path = os.path.join(output_dir, 'sequence.fasta')
        with open(seq_path, 'w') as f:
            f.write(f">Chromatogram_sequence\n{data['sequence']}\n")
        click.echo(f"Sequence saved to: {seq_path}")
        
    except Exception as e:
        click.echo(f"Error analyzing chromatogram: {e}", err=True)
        raise click.Abort()


@cli.command()
@click.argument('chromatogram_file', type=click.Path(exists=True))
@click.argument('reference_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', type=click.Path(), default='output',
              help='Output directory for results')
@click.option('--local-align', is_flag=True, 
              help='Perform local alignment instead of global')
def align(chromatogram_file: str, reference_file: str, output_dir: str, local_align: bool):
    """
    Align chromatogram sequence to a reference genome.
    
    CHROMATOGRAM_FILE: Path to the chromatogram file
    REFERENCE_FILE: Path to reference genome file (FASTA format)
    """
    click.echo(f"Aligning {chromatogram_file} to reference {reference_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read chromatogram
        reader = ChromatogramReader(chromatogram_file)
        query_sequence = reader.sequence
        
        # Read reference
        with open(reference_file, 'r') as f:
            lines = f.readlines()
            reference_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        
        # Perform alignment
        aligner = SequenceAligner(reference_sequence)
        alignment_result = aligner.align(query_sequence, local=local_align)
        
        # Report results
        click.echo(f"Alignment score: {alignment_result['score']:.2f}")
        click.echo(f"Sequence identity: {alignment_result['identity']:.2f}%")
        click.echo(f"Number of mismatches: {len(alignment_result['mismatches'])}")
        click.echo(f"Number of gaps: {alignment_result['gaps']}")
        
        # Save alignment results
        alignment_path = os.path.join(output_dir, 'alignment.txt')
        with open(alignment_path, 'w') as f:
            f.write("ALIGNMENT RESULTS\n")
            f.write("================\n\n")
            f.write(f"Score: {alignment_result['score']:.2f}\n")
            f.write(f"Identity: {alignment_result['identity']:.2f}%\n")
            f.write(f"Mismatches: {len(alignment_result['mismatches'])}\n")
            f.write(f"Gaps: {alignment_result['gaps']}\n\n")
            f.write("ALIGNMENT:\n")
            f.write(f"Reference: {alignment_result['aligned_reference']}\n")
            f.write(f"Query:     {alignment_result['aligned_query']}\n\n")
            
            if alignment_result['mismatches']:
                f.write("MISMATCHES:\n")
                for mm in alignment_result['mismatches'][:10]:  # Show first 10
                    f.write(f"Position {mm['position']}: {mm['reference']} -> {mm['query']}\n")
        
        click.echo(f"Alignment results saved to: {alignment_path}")
        
    except Exception as e:
        click.echo(f"Error performing alignment: {e}", err=True)
        raise click.Abort()


@cli.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.option('--reference-file', '-r', type=click.Path(),
              help='Reference genome file for alignment')
@click.option('--output-dir', '-o', type=click.Path(), default='batch_output',
              help='Output directory for batch results')
def batch(input_dir: str, reference_file: Optional[str], output_dir: str):
    """
    Process multiple chromatogram files in batch mode.
    
    INPUT_DIR: Directory containing chromatogram files
    """
    click.echo(f"Processing chromatogram files in: {input_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find chromatogram files
    input_path = Path(input_dir)
    chromatogram_files = list(input_path.glob('*.ab1')) + list(input_path.glob('*.abi')) + list(input_path.glob('*.scf'))
    
    if not chromatogram_files:
        click.echo("No chromatogram files found in input directory", err=True)
        return
    
    click.echo(f"Found {len(chromatogram_files)} chromatogram files")
    
    # Process each file
    for i, file_path in enumerate(chromatogram_files, 1):
        click.echo(f"Processing {i}/{len(chromatogram_files)}: {file_path.name}")
        
        try:
            # Create subdirectory for this file
            file_output_dir = os.path.join(output_dir, file_path.stem)
            os.makedirs(file_output_dir, exist_ok=True)
            
            # Analyze chromatogram
            reader = ChromatogramReader(str(file_path))
            data = reader.read()
            
            # Generate basic report
            report_path = os.path.join(file_output_dir, 'report.txt')
            with open(report_path, 'w') as f:
                f.write(f"CHROMATOGRAM ANALYSIS REPORT\n")
                f.write(f"============================\n\n")
                f.write(f"File: {file_path.name}\n")
                f.write(f"Sequence length: {len(data['sequence'])} bases\n")
                f.write(f"Average quality: {sum(data['quality_scores'])/len(data['quality_scores']):.2f}\n")
                f.write(f"Sequence:\n{data['sequence']}\n")
            
        except Exception as e:
            click.echo(f"Error processing {file_path.name}: {e}", err=True)
            continue
    
    click.echo(f"Batch processing complete. Results saved to: {output_dir}")


if __name__ == '__main__':
    cli()