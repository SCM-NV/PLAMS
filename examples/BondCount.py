#!/usr/bin/env amspython

from scm.plams import *
import matplotlib.pyplot as plt
import sys, os
import numpy as np

"""
Example to perform bond count based on a reactive MD simulation.

For each frame, plot the number of bond for every bond type in the simulation.

Bond type are based on rounded bond order (either integer or half-integer).

Usage:
$AMSBIN/amspython BondCount.py /path/to/ams.rkf
"""

def round_off(number,rtype='integer'):
    """
       Round off bond order to integer or half integer.
    """
    if rtype=='half_integer':
        number = round(number * 2) / 2
    if rtype=='integer':
        number = int(np.ceil(number)) 
    return number 

def bo2symbol(bo):
    """
       Convert bond order to bond symbol.
    """
    s = '---'
    if bo == 0.5: s = '--'
    if bo == 1.0: s = '-'
    if bo == 1.5: s = '=='
    if bo == 2.0: s = '='
    if bo == 2.5: s = '≡≡'
    if bo == 3.0: s = '≡'
    return s

def analyze(ams_rkf_path,round_type):
    """
        ams_rkf_path: str
            Path to an ams.rkf file from a reactive MD simulation

        Returns: dict (number of bond type)
            Each element of the dict are a list of length nframes containing 
            the occurence of every bond type at a given frame. 
    """

    if not os.path.exists(ams_rkf_path):
        raise FileNotFoundError(f"Couldn't find the file {ams_rkf_path}")

    # Load the rkf file
    job = AMSJob.load_external(ams_rkf_path)
    # Get number of frames
    ene_list = job.results.get_history_property('Energy')
    nframes = len(ene_list)+1

    # Loop over frames
    allbonds = {}
    for iframe in range(1,nframes):
        imol = job.results.get_history_molecule(iframe)
        bonds,btype = [],[]

        # Loop over atoms
        for i, atom in enumerate(imol.atoms, 1):

            # Loop over bonds to atom i
            for j, bond in enumerate(atom.bonds, 1):
                # Bond indices and symbols 
                b_atoms = [bond.atom1,bond.atom2]
                b_symbols = [bond.atom1.symbol,bond.atom2.symbol]
                sorted_symbols = sorted(b_symbols)

                # If bond has not yet been counted
                if b_atoms not in bonds:
                    # Compute BO and add bond to list
                    sorted_symbols.append(round_off(bond.order,round_type))
                    btype.append(sorted_symbols)
                    bonds.append(b_atoms)
                    bonds.append([b_atoms[1],b_atoms[0]])

        # Loop over unique bonds
        btype_set = set(tuple(row) for row in btype)
        for ibtype in btype_set:
            # Count how many bond of given type
            n_ibond = btype.count(list(ibtype))
            b_tag = ibtype[0]+bo2symbol(ibtype[2])+ibtype[1]
            b_label = '{:s} (BO {:2.1f})'.format(b_tag,ibtype[-1])
            # If new bond at frame iframe
            #    then add zeros for all frames before
            if b_label not in allbonds:
                allbonds[b_label] = [0.0]*(nframes-1)
                allbonds[b_label][iframe-1] = n_ibond
            # Else add bond
            else:
                allbonds[b_label][iframe-1] = n_ibond

        # Counter
        print('Frame {:d}/{:d}'.format(iframe,nframes))

    return allbonds

def plot_results(allbonds):
    """
        Plot all bonds vs. frame
    """
    for ibond in allbonds:
        plt.plot(allbonds[ibond],label=ibond)
    plt.legend()
    plt.xlabel('Frame')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig('bonds_analysis.pdf')
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: $AMSBIN/amspython bond.py /path/to/ams.rkf")
        exit(1)

    ams_rkf_path = sys.argv[1]

    # integer or half-integer
    round_type = 'integer'

    try:
        allbonds = analyze(ams_rkf_path,round_type)
        plot_results(allbonds)
    
    except Exception as e:
        print(e)
        exit(1)

if __name__ == '__main__':
    main()
