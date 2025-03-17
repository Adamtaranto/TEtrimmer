from os import path
from typing import Optional

import click
from derip2.derip import DeRIP


def get_derip_consensus(
    input_file: str,
    output_file: str,
    consensus_name: str = 'derip_consensus',
    maxSNPnoise: float = 0.5,
    minRIPlike: float = 0.1,
    maxGaps: float = 0.7,
    reaminate: bool = False,
    fillindex: Optional[int] = None,
    fillmaxgc: bool = False,
):
    """
    Generate a deRIPed consensus sequence from an alignment file.

    This function processes a DNA alignment file to identify and correct Repeat-Induced
    Point (RIP) mutations, then writes the corrected consensus sequence to an output file.
    It also prints a summary of identified RIP mutations and alignment information.

    Parameters
    ----------
    input_file : str
        Path to the input alignment file in FASTA format.
    output_file : str
        Path where the consensus sequence will be written.
    consensus_name : str, optional
        Name for the consensus sequence in the output file (default: 'derip_consensus').
    maxSNPnoise : float, optional
        Maximum proportion of conflicting SNPs permitted before excluding a column
        from RIP/deamination assessment (default: 0.5).
    minRIPlike : float, optional
        Minimum proportion of deamination events in RIP context required for a
        column to be deRIPed in final sequence (default: 0.1).
    maxGaps : float, optional
        Maximum proportion of gaps in a column before considering it a gap
        in the consensus sequence (default: 0.7).
    reaminate : bool, optional
        Whether to correct all deamination events independent of RIP context (default: False).
    fillindex : int, optional
        Index of row to use for filling uncorrected positions (default: None).
    fillmaxgc : bool, optional
        Whether to use the sequence with highest GC content for filling if
        no row index is specified (default: False).

    Returns
    -------
    None
        This function writes to a file and prints to stdout but doesn't return a value.

    Raises
    ------
    FileNotFoundError
        If the specified input file does not exist.

    Notes
    -----
    The function prints summary information to standard output, including:
    - Number of columns repaired
    - RIP mutation summary by sequence
    - Visualization of the masked alignment
    - Gapped consensus sequence
    """
    if path.isfile(input_file):
        derip_object = DeRIP(
            input_file,
            maxSNPnoise=maxSNPnoise,
            minRIPlike=minRIPlike,
            maxGaps=maxGaps,
            reaminate=reaminate,
            fillindex=fillindex,
            fillmaxgc=fillmaxgc,
        )
        # Calculate RIP mutations
        derip_object.calculate_rip(label=consensus_name)

        # Access corrected positions
        print(
            f'DeRIP2 found {len(derip_object.corrected_positions)} columns to be repaired.\n'
        )

        # Print RIP summary
        print('RIP summary by row:')
        derip_object.print_rip_summary()

        # Print masked alignment
        print('\nMutation masked alignment:\n', derip_object.masked_alignment)
        print(f'{str(derip_object.gapped_consensus.seq)} {consensus_name}\n')

        # Opt1: Write output consensus to file
        derip_object.write_consensus(output_file, consensus_id=consensus_name)

        # Opt2: Write original alignment with appened deRIP'd sequence to output file
        # derip_object.write_alignment(output_file, append_consensus=True, mask_rip=False)

    else:
        raise FileNotFoundError(f"The file '{input_file}' does not exist.")


@click.command()
@click.option(
    '--input_file',
    '-i',
    required=True,
    type=str,
    help='Multiple sequence alignment FASTA file path',
)
@click.option('--output_file', '-o', required=True, type=str, help='Output file')
@click.option(
    '--consensus_name',
    '-n',
    type=str,
    default='derip_consensus',
    help='Name of the consensus sequence (default: derip_consensus)',
)
@click.option(
    '--maxSNPnoise',
    type=float,
    default=0.5,
    help='Maximum proportion of conflicting SNPs permitted before excluding column from RIP/deamination assessment (default: 0.5)',
)
@click.option(
    '--minRIPlike',
    type=float,
    default=0.1,
    help='Minimum proportion of deamination events in RIP context required for column to be deRIPd in final sequence (default: 0.1)',
)
@click.option(
    '--maxGaps',
    type=float,
    default=0.7,
    help='Maximum proportion of gaps in a column before considering it a gap in consensus (default: 0.7)',
)
@click.option(
    '--reaminate',
    is_flag=True,
    help='Correct all deamination events independent of RIP context',
)
@click.option(
    '--fillindex',
    type=int,
    help='Index of row to use for filling uncorrected positions',
)
@click.option(
    '--fillmaxgc',
    is_flag=True,
    help='Use sequence with highest GC content for filling if no row index is specified',
)
def derip_click(
    input_file,
    output_file,
    consensus_name,
    maxsnpnoise,
    minriplike,
    maxgaps,
    reaminate,
    fillindex,
    fillmaxgc,
):
    """
    Command line interface for the deRIP consensus generation tool.

    This function serves as the entry point for the command line interface,
    processing arguments from Click decorators and passing them to the
    get_derip_consensus function.

    Parameters
    ----------
    input_file : str
        Path to the input alignment file in FASTA format.
    output_file : str
        Path where the consensus sequence will be written.
    consensus_name : str
        Name for the consensus sequence in the output file.
    maxsnpnoise : float
        Maximum proportion of conflicting SNPs permitted before excluding a column
        from RIP/deamination assessment.
    minriplike : float
        Minimum proportion of deamination events in RIP context required for a
        column to be deRIPed in final sequence.
    maxgaps : float
        Maximum proportion of gaps in a column before considering it a gap
        in the consensus sequence.
    reaminate : bool
        Whether to correct all deamination events independent of RIP context.
    fillindex : int or None
        Index of row to use for filling uncorrected positions.
    fillmaxgc : bool
        Whether to use the sequence with highest GC content for filling if
        no row index is specified.

    Returns
    -------
    None
        This function calls get_derip_consensus which writes output to a file
        and prints information to stdout.

    Notes
    -----
    This function is intended to be used with Click as a command-line entry point
    and should not typically be called directly in code.
    """
    get_derip_consensus(
        input_file=input_file,
        output_file=output_file,
        consensus_name=consensus_name,
        maxSNPnoise=maxsnpnoise,
        minRIPlike=minriplike,
        maxGaps=maxgaps,
        reaminate=reaminate,
        fillindex=fillindex,
        fillmaxgc=fillmaxgc,
    )


if __name__ == '__main__':
    derip_click()
