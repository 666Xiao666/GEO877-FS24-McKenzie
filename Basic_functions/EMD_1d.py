def first_nonzero(vec):
    """
    Finds the index of the first non-zero element in a vector.

    Args:
    vec (list): The input vector.

    Returns:
    int: The index of the first non-zero element, or -1 if no non-zero element is found.
    """
    dim = len(vec)
    for i in range(dim):
        if vec[i] > 0.0:
            return i
    return -1  # no empty cells found

def move_dirt(dirt, di, holes, hi):
    """
    Moves dirt from one location to another, calculating the work done.

    Args:
    dirt (list): The amount of dirt at each location.
    di (int): Index of the source location.
    holes (list): The amount of available space at each destination.
    hi (int): Index of the destination location.

    Returns:
    float: The work done to move the dirt.
    """
    # Move as much dirt as possible from di to hi
    if dirt[di] <= holes[hi]:   # use all dirt
        flow = dirt[di]
        dirt[di] = 0.0            # all dirt got moved
        holes[hi] -= flow         # less to fill now
    elif dirt[di] > holes[hi]:  # use just part of dirt
        flow = holes[hi]          # fill remainder of hole
        dirt[di] -= flow          # less dirt left
        holes[hi] = 0.0           # hole is filled
    dist = abs(di - hi)  # Calculate the distance
    return flow * dist    # Calculate the work done

def my_wasserstein(p, q):
    """
    Computes the Wasserstein distance between two distributions.

    Args:
    p (list): The first distribution.
    q (list): The second distribution.

    Returns:
    float: The Wasserstein distance between the distributions.
    """
    # Create copies of p and q
    dirt = p[:]
    holes = q[:]
    tot_work = 0.0

    # Iterate until all dirt is moved or all holes are filled
    while True:  
        from_idx = first_nonzero(dirt)
        to_idx = first_nonzero(holes)
        if from_idx == -1 or to_idx == -1:
            break
        # Move dirt from the source to the destination and calculate work done
        work = move_dirt(dirt, from_idx, holes, to_idx)
        tot_work += work
    return tot_work  

def main():
    print("\nBegin Wasserstein distance demo ")

    # Define distributions 
    P = [0.6, 0.1, 0.1, 0.1, 0.1]
    Q1 = [0.1, 0.1, 0.6, 0.1, 0.1]
    Q2 = [0.1, 0.1, 0.1, 0.1, 0.6]

    # Compute Wasserstein distances
    wass_p_q1 = my_wasserstein(P, Q1)
    wass_p_q2 = my_wasserstein(P, Q2)

    # Print results
    print("\nWasserstein distances: ")
    print("P to Q1 : %0.4f " % wass_p_q1)
    print("P to Q2 : %0.4f " % wass_p_q2)

    print("\nEnd demo ")

if __name__ == "__main__":
    main()
