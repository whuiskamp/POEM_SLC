def flood_fill(i, j, wet, omask, xcyclic=True, tripolar=True):
    #
    # An iterative (stack based) flood fill algorithm. This is a modified
    # version of the algorithm found in NCAR's topo edit notebook.
    # Link: https://github.com/NCAR/tx2_3/blob/main/topography/MaskEdit_tx2_3v2b.ipynb
    # 
    # The flood fill starts at [j,i] and treats any positive value of "omask" as
    # passable. Zero and negative values block flooding.

    # xcyclic = True allows cyclic behavior in the last index. (default)
    # tripolar = True allows a fold across the top-most edge. (default)
    #
    # Returns an array of 0's (not wet) and positive integers (defined by 'wet' - wet).
    #
    iso_mask = 0*omask
    (nj,ni) = iso_mask.shape
    stack = set()
    stack.add( (j,i) )
    while stack:
        (j,i) = stack.pop()
        if iso_mask[j,i] or omask[j,i] <= 0: continue
        iso_mask[j,i] = wet
        if i>0: stack.add( (j,i-1) )
        elif xcyclic: stack.add( (j,ni-1) )
        if i<ni-1: stack.add( (j,i+1) )
        elif xcyclic: stack.add( (j,0) )
        if j>0: stack.add( (j-1,i) )
        if j<nj-1: stack.add( (j+1,i) )
        elif tripolar: stack.add( (j,ni-1-i) ) # Tri-polar fold
    return iso_mask