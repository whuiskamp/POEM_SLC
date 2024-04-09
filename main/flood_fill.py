def flood_fill(col,row,ocean,filled,iso_mask):
    # Recursive function to fill in the ocean and, as a result
    # identify isolated cells created during the wetting/drying 
    # in the previous steps.
    # row,col  = Initial row/col values which we are certain are open ocean
    # ocean    = Value of open ocean in the o_mask, in this case, 1
    # filled   = Value of filled ocean, in this case 0
    # iso_mask = The field we want to run on. Should be a (deep)copy
    #              of the new ocean mask.  
    # 
    # Check we're inside the model domain
    if col < 0 or col >= len(iso_mask[1]) or row < 0 or row >= len(iso_mask[0]):
        return
    # Is the current cell ocean? If not, return
    if iso_mask[row,col] != ocean:
        return
    # If so, fill
    iso_mask[row,col] = filled
    # fourthly, attempt to fill the neighboring positions
    flood_fill(col+1,row,ocean,filled,iso_mask)
    flood_fill(col-1,row,ocean,filled,iso_mask)
    flood_fill(col,row+1,ocean,filled,iso_mask)
    flood_fill(col,row-1,ocean,filled,iso_mask)
