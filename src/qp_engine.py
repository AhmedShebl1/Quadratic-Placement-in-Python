import numpy as np
import sys
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

MIN_Y = MIN_X = 0
MAX_Y = MAX_X = 100
CHIP_DIMENSIONS = ((MIN_X, MAX_X), (MIN_Y, MAX_Y))

global count
count = 1

def parse_input_file(file_path):
    """
    File format :

    NumberofGates G and NumberofNets N
    Gate# NumNetsConnected NetNumber...
    ...
    Number of pads P connected to nets
    PinID NetNumberConnectedTo PinX PinY
    ...

    Example :    

    18 20
    1 2 20 19
    ...
    6
    1 12 25 0
    ...

    """

    # Initialize variables 
    num_gates = num_pads = num_nets = 0

    # Includes the nets connected to each gate
    gate_net_map = {}

    # Includes the pad's coordinates and connected net 
    pad_info = {}

    # Includes info about each net and its connected gates and pads
    nets = {}

    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        num_gates, num_nets = map(int, first_line.split())

        for _ in range(num_gates):
            line = file.readline().strip()
            parts = list(map(int, line.split()))
            gate_id = parts[0]
            nets_connected = parts[2:]
            gate_net_map[gate_id] = nets_connected

        num_pads = int(file.readline().strip())

        for _ in range(num_pads):
            line = file.readline().strip()
            parts = list(map(int, line.split()))
            pin_id = parts[0]
            net_connected = parts[1]
            pin_x = parts[2]
            pin_y = parts[3]
            pad_info[pin_id] = {
                'net': net_connected,
                'x': pin_x,
                'y': pin_y
            }

    # Get the gates and pads connected to each net
    for net_id in range(1, num_nets + 1):
        nets[net_id] = {
            'gates': [],
            'pads': []
        }
        for gate_id, connected_nets in gate_net_map.items():
            if net_id in connected_nets:
                nets[net_id]['gates'].append(gate_id)
        for pad_id, pad_details in pad_info.items():
            if pad_details['net'] == net_id:
                nets[net_id]['pads'].append(pad_id)

    #Initialize gate positions
    initial_gate_pos = {}
    
    # Initial random placement of gates
    for g_id in range(1 , num_gates + 1):
        initial_gate_pos[g_id] = (np.random.randint(0, MAX_X), np.random.randint(0, MAX_Y))


    return initial_gate_pos, gate_net_map, pad_info, nets

def gen_output_file (gate_pos, file_path):
    """
    Generates an output file with the final gate positions.
    The format is:
    GateID X Y
    positions are written to 8 decimal places
    """
    with open(file_path, 'w') as file:
        for gate_id, (x, y) in gate_pos.items():
            file.write(f"{gate_id} {x:.8f} {y:.8f}\n")
    
    print(f"Output written to {file_path}")

def recursive_place(gate_ids, current_dim, local_gate_pos, gate_net_map, pad_info, nets, depth_remaining, slice_direction = True):
    """
    This helper function places gates recursively for the divided subreigons based on the QP placement algorithm function
    The placer slices the area into vertical and horizontal regions until the assigned depth is reached.
    """

    # declare gae_pos global because we mutate it below
    global gate_pos, count

    # Debug: Print current iteration
    print("Iteration:", count)
    count += 1

    # Guard against empty gate list
    if not gate_ids:
        return
    
    # Perform placement for the current region
    new_gate_pos = qp_place(gate_ids, local_gate_pos, gate_net_map, pad_info, nets)
    local_gate_pos.update(new_gate_pos)

    if depth_remaining == 0:
        # Placement is complete
        # Edit the global gate_pos with final positions for this region's gates
        for gate_id in gate_ids:
            gate_pos[gate_id] = local_gate_pos[gate_id]
        return

    # Unpack the dimensions
    (min_x, max_x), (min_y, max_y) = current_dim

    # Initialize variables for slicing (left/bottom) for the first set and (right/top) the second set of gates 
    gate_ids_1 = gate_ids_2 = []
    dim_1 = dim_2 = ()
    pad_info_1 = pad_info_2 = {}
    nets_1 = nets_2 = {}
    local_gate_pos_1 = local_gate_pos_2 = {}
    gate_net_map_1 = gate_net_map_2 = {}
    mid_pt = 0

    # Slice the area in the specified direction
    if slice_direction:
        # Vertical slice
        # Sort gates based on x-coordinates and split
        gate_ids_1, gate_ids_2 = slice_chip(gate_ids, local_gate_pos, axis=0)

        # Update current_dim for left and right halves
        mid_pt = (min_x + max_x) / 2
        dim_1 = ((min_x, mid_pt), (min_y, max_y))
        dim_2 = ((mid_pt, max_x), (min_y, max_y))
        pass
    else:
        # Horizontal slice
        # Sort gates based on y-coordinates and split
        gate_ids_1, gate_ids_2 = slice_chip(gate_ids, local_gate_pos, axis=1)

        # Update current_dim for upper and lower halves
        mid_pt = (min_y + max_y) / 2
        dim_1 = ((min_x, max_x), (mid_pt, max_y))
        dim_2 = ((min_x, max_x), (min_y, mid_pt))

        pass
    
    # Propagate connected gates and pads outside the current region to the edges
    pad_info_1, gate_net_map_1, local_gate_pos_1 = propagate_to_edges(gate_ids_1, local_gate_pos, gate_net_map, pad_info, dim_1)
    pad_info_2, gate_net_map_2, local_gate_pos_2 = propagate_to_edges(gate_ids_2, local_gate_pos, gate_net_map, pad_info, dim_2)

    # Compute the net connections with the newly generated psuedo-pads
    nets_1 = compute_nets(gate_ids_1, gate_net_map_1, pad_info_1)
    nets_2 = compute_nets(gate_ids_2, gate_net_map_2, pad_info_2)

    # Recur for the next depth level on both halfs and alternate the slicing direction
    recursive_place(gate_ids_1, dim_1, local_gate_pos_1, gate_net_map_1, pad_info_1, nets_1, depth_remaining - 1, not slice_direction)
    recursive_place(gate_ids_2, dim_2, local_gate_pos_2, gate_net_map_2, pad_info_2, nets_2, depth_remaining - 1, not slice_direction)

    return        

def propagate_to_edges(gate_ids, local_gate_pos, gate_net_map, pad_info, current_dim):
    """
    Propagates the positions of gates and pads outside the current subregion
    to the edges of the subregion.
    """

    # Unpack region boundaries
    (min_x, max_x), (min_y, max_y) = current_dim

    # Collect all nets connected to gates inside this region
    connected_nets = set()
    for gate_id in gate_ids:
        for net in gate_net_map.get(gate_id, []):
            connected_nets.add(net)

    # Build a unified list of all external objects (gates + pads)
    external_objects = []

    # External gates: those NOT in the subregion
    for gate_id in gate_net_map:
        if gate_id not in gate_ids:
            external_objects.append(("gate", gate_id))

    # All pads are added, but we will **not project pads inside the region**
    for pad_id in pad_info:
        external_objects.append(("pad", pad_id))

    # For each external object, project it onto the boundary 
    # Structure : projected_positions[(kind, id)] = (x_proj, y_proj)
    projected_positions = {}

    for kind, obj_id in external_objects:

        # Get this object's nets
        if kind == "gate":
            nets_list = gate_net_map[obj_id]
            x_prev, y_prev = local_gate_pos[obj_id]

            # Skip if this external gate does not touch region
            if set(nets_list).isdisjoint(connected_nets):
                continue

        else:  # pad
            nets_list = [pad_info[obj_id]['net']]
            x_prev = pad_info[obj_id]['x']
            y_prev = pad_info[obj_id]['y']

            # Pads INSIDE the region must NOT be projected
            if min_x <= x_prev <= max_x and min_y <= y_prev <= max_y:
                continue

            # Skip pads not connected to this region
            if pad_info[obj_id]['net'] not in connected_nets:
                continue

        # Determine which boundary point to collapse to

        # Project strictly onto region boundary
        if x_prev < min_x:
            x_new = min_x
        elif x_prev > max_x:
            x_new = max_x
        else:
            x_new = x_prev

        if y_prev < min_y:
            y_new = min_y
        elif y_prev > max_y:
            y_new = max_y
        else:
            y_new = y_prev


        # Store projected position
        projected_positions[(kind, obj_id)] = (x_new, y_new)

    # Return the projected positions as pseudo-pads
    # Keep the pads inside the region unchanged
    extended_pad_info = {}
    next_pad_id = 1

    # Keep pads inside the region ONLY if their net touches an internal gate
    for pad_id, info in pad_info.items():
        net = info['net']
        x, y = info['x'], info['y']

        if not (min_x <= x <= max_x and min_y <= y <= max_y):
            continue

        # Skip pads whose net is not connected to any internal gate
        if net not in connected_nets:
            continue

        extended_pad_info[next_pad_id] = {
            'net': net,
            'x': x,
            'y': y
        }
        next_pad_id += 1


    # Add projected objects as pseudo-pads
    for (kind, obj_id), (x_proj, y_proj) in projected_positions.items():

        if kind == "pad":
            # Update existing pad position
            extended_pad_info[next_pad_id] = {
                'net': pad_info[obj_id]['net'],
                'x': x_proj,
                'y': y_proj
            }
            next_pad_id += 1

        else:
            # Create a new pseudo-pad for the external gate 
            # For each net connected to this gate, create a pseudo-pad
            for net in gate_net_map[obj_id]:
                # Skip nets not connected to internal gates
                if net not in connected_nets:
                    continue
                extended_pad_info[next_pad_id] = {
                    'net': net,
                    'x': x_proj,
                    'y': y_proj
                }
                next_pad_id += 1

    # Return the new gate_net_map and local_gate_pos with only the internal gates
    new_gate_net_map = {}
    new_local_gate_pos = {}
    for gate_id in gate_ids:
        new_gate_net_map[gate_id] = gate_net_map[gate_id]
        new_local_gate_pos[gate_id] = local_gate_pos[gate_id]

    # Return these projected boundary anchors which will be treated as pseudo-pads for the next QP solve.
    return extended_pad_info, new_gate_net_map, new_local_gate_pos

def compute_nets(gate_ids, gate_net_map, pad_info):
    """
    Computes the nets relevant to this subregion.
    Only nets that have at least one internal gate are included.
    Pads are included only if their net is in this set.
    """
    nets = {}

    # Collect nets that internal gates use
    internal_nets = set()
    for gate_id in gate_ids:
        for net in gate_net_map.get(gate_id, []):
            internal_nets.add(net)

    # Initialize net structures with empty lists
    for net in internal_nets:
        nets[net] = {'gates': [], 'pads': []}

    # Fill gate lists
    for gate_id in gate_ids:
        for net in gate_net_map.get(gate_id, []):
            nets[net]['gates'].append(gate_id)

    # Add only pads for nets used by internal gates
    for pad_id, pad_details in pad_info.items():
        net = pad_details['net']
        if net in internal_nets:
            nets[net]['pads'].append(pad_id)

    return nets

def slice_chip(gate_ids, local_gate_pos, axis=0):
    """
    Sorts the gates based on their positions along the specified axis.
    axis=0 for x-axis, axis=1 for y-axis
    """

    sorted_gates = sorted(gate_ids, key=lambda gid: local_gate_pos[gid][axis])
    mid = len(sorted_gates) // 2

    return sorted_gates[:mid], sorted_gates[mid:] 

def get_connectivity_matrix(gate_ids, gate_net_map, nets):
    """
    Constructs the connectivity matrix for the given gates and pads.
    """
    mat = np.zeros((len(gate_ids), len(gate_ids)))

    for idx1, i in enumerate(gate_ids):
        for idx2, j in enumerate(gate_ids):
            if i != j:
                # Compute sum of clique weights for all nets shared by i and j.
                # Each net of size k contributes weight = 1/(k-1) to every pair on that net.
                shared = set(gate_net_map.get(i, [])).intersection(set(gate_net_map.get(j, [])))
                w = 0.0
                if shared:
                    for net in shared:
                        # k = total members of the net (gates + pads)
                        k = len(nets[net]['gates']) + len(nets[net]['pads'])
                        if k <= 1:
                            # degenerate net: skip to avoid div-by-zero
                            continue
                        w += 1.0 / (k - 1)
                mat[idx1, idx2] = w
            else :
                mat[idx1][idx2] = 0  # No self-connections
    
    return mat

def compute_A_matrix(connectivity_matrix, pad_info, gate_ids, gate_net_map, nets):
    """
    Computes the A matrix for the QP formulation.
    """

    # Determine pad connections for each gate
    pad_connections = np.zeros(connectivity_matrix.shape[0])
    for pad_id in pad_info:
        net = pad_info[pad_id]['net']
        # compute k based on gates+pads on this net
        k = len(nets[net]['gates']) + len(nets[net]['pads'])
        if k <= 1:
            continue
        w = 1.0 / (k - 1)

        # For each gate in this subproblem,
        # if this gate is connected to the same net as the pad,
        # add the pad’s pull weight to that gate’s diagonal entry
        for idx, gate_id in enumerate(gate_ids):
            if net in gate_net_map.get(gate_id, []):
                pad_connections[idx] += w

    # Diagonal elements are the sum of connections for each gate and pad connections
    degree_matrix = np.diag(np.sum(connectivity_matrix, axis=1) + pad_connections) 

    # Non-diagonal elements are negative connections
    A = degree_matrix - connectivity_matrix

    return A

def compute_b_vector(gate_ids, gate_net_map, pad_info, nets):
    """
    Computes the b vector for the QP formulation.
    """
    b_x = np.zeros(len(gate_ids))
    b_y = np.zeros(len(gate_ids))

    # For each gate, sum the positions of connected pads
    for idx, gate_id in enumerate(gate_ids):
        nets_list = gate_net_map.get(gate_id, [])
        for net in nets_list:
            # weight pads by 1/(k-1) where k is total members (gates + pads)
            k = len(nets[net]['gates']) + len(nets[net]['pads'])
            if k <= 1:
                continue
            w = 1.0 / (k - 1)
            for pad_id in pad_info:
                if pad_info[pad_id]['net'] == net:
                    b_x[idx] += pad_info[pad_id]['x'] * w
                    b_y[idx] += pad_info[pad_id]['y'] * w
    
    return b_x, b_y

def qp_place(gate_ids, local_gate_pos, gate_net_map, pad_info, nets):
    """
    Solve QP for the given subset of gate_ids and update local_gate_pos.
    """

    # Handdle empty gate list
    if not gate_ids:
        return {}

    # Get the connectivity matrix
    connectivity_matrix = get_connectivity_matrix(gate_ids, gate_net_map, nets)

    # Derive the A matrix from the connectivity matrix
    A = compute_A_matrix(connectivity_matrix, pad_info, gate_ids, gate_net_map, nets)

    # Get the b vectors for x and y coordinates
    b_x, b_y = compute_b_vector(gate_ids, gate_net_map, pad_info, nets)

    # Solve for optimal positions using sparse solver
    A_sparse = coo_matrix(A).tocsr()
    x_positions = spsolve(A_sparse, b_x)
    y_positions = spsolve(A_sparse, b_y)
    
    # Update local_gate_pos with new positions
    for idx, gate_id in enumerate(gate_ids):
        local_gate_pos[gate_id] = (x_positions[idx], y_positions[idx])

    return local_gate_pos

def qp_engine(file_path,  depth=1):
    """
    Main function to run the QP placement algorithm. Calls parsing, recursive placement, and output generation.
    """
 
    # Parse the input file to get gate positions, net connections, and pad info
    global gate_pos
    gate_net_map, pad_info, nets = {}, {}, {}
    gate_pos, gate_net_map, pad_info, nets = parse_input_file(file_path)

    # Start the recursive placement process
    recursive_place(list(gate_pos.keys()), CHIP_DIMENSIONS, gate_pos, gate_net_map, pad_info, nets, depth, True)

    # Generate the output file with final gate positions
    gen_output_file(gate_pos, './outputs_zero_initial/' + file_path.split('/')[-1])

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("Usage: python qp_engine.py <input_file_path> <depth>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    target_depth = int(sys.argv[2])
    qp_engine(input_file_path, target_depth)
