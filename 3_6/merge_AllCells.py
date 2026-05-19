#!/usr/bin/env python3
"""Merge per-block AllCells files into a single Tecplot FEPOINT BRICK file.

Usage:
    python merge_AllCells.py <AllCells_XXXXXXXX>   # single directory
    python merge_AllCells.py -a                    # all AllCells_* in cwd
"""
import sys
import os
import glob


def merge(dirpath):
    dirpath = dirpath.rstrip('/')
    donefile = os.path.join(dirpath, '.merged')
    if os.path.exists(donefile):
        print(f'Already processed: {dirpath} (remove {donefile} to reprocess)')
        return

    # --- Parse info.dat ---
    info = {}
    variables_line = None
    with open(os.path.join(dirpath, 'info.dat')) as f:
        for line in f:
            line = line.strip()
            if line.startswith('VARIABLES'):
                variables_line = line
            elif '=' in line:
                key, val = line.split('=', 1)
                info[key.strip()] = val.strip()

    nBlocks = int(info['nBlocks'])
    ni      = int(info['ni'])
    nj      = int(info['nj'])
    nk      = int(info['nk'])

    nCorners         = nBlocks * (ni+1) * (nj+1) * (nk+1)
    nCells           = nBlocks * ni * nj * nk
    nCornersPerBlock = (ni+1) * (nj+1) * (nk+1)

    # Strides for 1-indexed corner (i,j,k) within a block:
    # index = (k-1)*(ni+1)*(nj+1) + (j-1)*(ni+1) + i
    stride_j = ni + 1
    stride_k = (ni+1) * (nj+1)

    outpath = os.path.join(dirpath, 'merged.dat')

    # Sort block files by their numeric block index
    block_files = sorted(
        [f for f in os.listdir(dirpath) if f.startswith('Block_') and f.endswith('.dat')],
        key=lambda f: int(f[6:14])
    )

    if len(block_files) != nBlocks:
        print(f"Warning: expected {nBlocks} block files, found {len(block_files)}")

    with open(outpath, 'w') as out:

        # --- Header ---
        out.write('TITLE = "AllCells"\n')
        out.write(variables_line + '\n')
        out.write(f'ZONE N={nCorners}, E={nCells}, F=FEPOINT, ET=BRICK\n')

        # --- Data: concatenate block files in block-index order ---
        for bfile in block_files:
            with open(os.path.join(dirpath, bfile)) as bf:
                out.write(bf.read())

        # --- Connectivity ---
        # Block b (0-indexed) occupies corners with 1-based offset b*nCornersPerBlock.
        # For cell (i,j,k) inside that block the 8 corners in Tecplot BRICK order are:
        #   bottom face: c1 c2 c3 c4  (k layer)
        #   top face:    c5 c6 c7 c8  (k+1 layer, = bottom + stride_k)
        for b in range(nBlocks):
            offset = b * nCornersPerBlock
            for k in range(1, nk+1):
                for j in range(1, nj+1):
                    for i in range(1, ni+1):
                        c1 = offset + (k-1)*stride_k + (j-1)*stride_j + i
                        c2 = c1 + 1
                        c3 = c1 + stride_j + 1
                        c4 = c1 + stride_j
                        c5 = c1 + stride_k
                        c6 = c2 + stride_k
                        c7 = c3 + stride_k
                        c8 = c4 + stride_k
                        out.write(f'{c1} {c2} {c3} {c4} {c5} {c6} {c7} {c8}\n')

    for bfile in block_files:
        os.remove(os.path.join(dirpath, bfile))
    open(donefile, 'w').close()
    print(f'Written: {outpath}  (N={nCorners}, E={nCells})')


if len(sys.argv) < 2:
    print("Usage: merge_AllCells.py <AllCells_XXXXXXXX dir>")
    print("       merge_AllCells.py -a   # process all AllCells_* in cwd")
    sys.exit(1)

if sys.argv[1] == '-a':
    dirs = sorted(glob.glob('AllCells_*'))
    if not dirs:
        print('No AllCells_* directories found.')
        sys.exit(1)
    for d in dirs:
        if os.path.isdir(d):
            merge(d)
else:
    merge(sys.argv[1])
