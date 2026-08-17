import os
import gmsh

#  =============================================================================
#    This python script relies on GMSH's python API to easily and more quickly
#    process and record mesh information upon which Lionbolt relies.
#    
#    Following processing, the data is written to text files formatted for
#    ease of reading within Lionbolt.
#  =============================================================================

# This script is certainly going to be removed in the near future,
# with Lionbolt relying on Fortran-only code to perform this task.
# There are many issues with it, including debugging difficulties,
# creation of temporary files, etc.

def parser ():
    import argparse
    
    # Original intention was for this to be a user-friendly script... not at all necessary anymore...
    p = argparse.ArgumentParser()
    p.add_argument('-d', action='store_true', help='Process this as debug')
    p.add_argument('-b', type=str, help='The base directory of Lionbolt.')
    p.add_argument('-f', type=str, help='The name of the mesh file to process.')
    p.add_argument('-s', type=str, help='Scratch folder name.')
    
    args = p.parse_args()
    
    return args.d, args.b, args.f, args.s

# get_faces function. This enforces GMSH's node ordering system for local faces
def get_faces(nodes, element_type, sort=False):
    def create_face(node_list):
        return tuple(sorted(node_list) if sort else node_list)
    
    if element_type == 4:    # 4-node tetrahedron
        return [
            create_face([nodes[0], nodes[1], nodes[2]]),
            create_face([nodes[0], nodes[1], nodes[3]]),
            create_face([nodes[0], nodes[2], nodes[3]]),
            create_face([nodes[1], nodes[2], nodes[3]])
        ]
    elif element_type == 5:  # 8-node hexahedron
        return [
            create_face([nodes[0], nodes[1], nodes[2], nodes[3]]),
            create_face([nodes[4], nodes[5], nodes[6], nodes[7]]),
            create_face([nodes[0], nodes[1], nodes[5], nodes[4]]),
            create_face([nodes[2], nodes[3], nodes[7], nodes[6]]),
            create_face([nodes[0], nodes[3], nodes[7], nodes[4]]),
            create_face([nodes[1], nodes[2], nodes[6], nodes[5]])
        ]
    elif element_type == 6:  # 6-node prism
        return [
            create_face([nodes[0], nodes[1], nodes[2]]),
            create_face([nodes[3], nodes[4], nodes[5]]),
            create_face([nodes[0], nodes[1], nodes[4], nodes[3]]),
            create_face([nodes[1], nodes[2], nodes[5], nodes[4]]),
            create_face([nodes[2], nodes[0], nodes[3], nodes[5]])
        ]
    elif element_type == 7:  # 5-node pyramid
        return [
            create_face([nodes[0], nodes[1], nodes[2], nodes[3]]),
            create_face([nodes[0], nodes[1], nodes[4]]),
            create_face([nodes[1], nodes[2], nodes[4]]),
            create_face([nodes[2], nodes[3], nodes[4]]),
            create_face([nodes[3], nodes[0], nodes[4]])
        ]
    else:
        return []

def raw_nodes ():
    # Get the nodes from GMSH then just write them in an array
    # Will I need node tags eventually?
    _, node_coords, _ = gmsh.model.mesh.getNodes()
    meshnodes = [(node_coords[i], node_coords[i+1], node_coords[i+2])
                 for i in range(0, len(node_coords), 3)]
    
    # Now simply write to raw_nodes.txt
    fname = os.path.join(bpath, f'{scratch}TAPE%0_raw_nodes.txt')
    with open(fname, 'w') as f:
        # FORMAT:
        # First column is node index, then the next 3 columns are xyz
        f.write(f'{len(meshnodes)}\n')
        for i in range(0, len(meshnodes)):
            x = meshnodes[i][0]
            y = meshnodes[i][1]
            z = meshnodes[i][2]
            f.write(f'{i+1} {x} {y} {z}\n')

def elements (): 
    # Get the 3D elements (defined by the nodes that make up the element)
    # Format them in a list of dictionaries, element_dict. Each entry corresponds to an element type
    # They're read in as entities, defined as per GMSH
    element_dict = []
    for dim, entity_tag in gmsh.model.getEntities(3):
        element_types, element_tags_blocks, node_tags_blocks = gmsh.model.mesh.getElements(dim, entity_tag)
        for element_type, tags, nodes in zip(element_types, element_tags_blocks, node_tags_blocks):
            NK = gmsh.model.mesh.getElementProperties(element_type)[3]
            block = {
                'type': element_type,
                'entity_tag': entity_tag,
                'dim': dim,
                'elements': []
            }
            for i in range(0, len(nodes), NK):
                eid = tags[i // NK]
                node_ids = nodes[i:i+NK]
                block['elements'].append((eid, node_ids))
            
            element_dict.append(block)
    
    # Now shift the 3D elements in case they're not being indexed starting from 1
    SHIFT = min(el for block in element_dict if block['dim'] == 3
                   for el, _ in block['elements']) - 1
    
    for block in element_dict:
        if block['dim'] == 3:
            block['elements'] = [(e - SHIFT, nodes) for e, nodes in block['elements']]
    
    # Determine how many 3D elements there are in the mesh
    NE = 0
    for block in element_dict:
        elements = block['elements']
        NE = NE + len(elements)
    
    # Now write the elements
    fname = os.path.join(bpath, f'{scratch}TAPE%0_elements.txt')
    with open(fname, 'w') as f:
        # FORMAT:
        # First, write the number of elements
        # Then, go through each element type and write all of the elements
        # Much like the raw GMSH file, there is a header for every element type
        # which gives the element type index as well as the number of elements
        # in this type. This is necessary for reading in Lionbolt.
        f.write(f'{NE}\n')
        for block in element_dict:
            etype = block['type']
            elements = block['elements']
            f.write(f'{etype} {len(elements)}\n')
            for element, nodes in elements:
                nodesstr = ' '.join(map(str, nodes))
                f.write(f'{element} {nodesstr}\n')

def define_global_faces ():
    # global_face_dict gives, when given a set of sorted nodes that are in a face,
    #   the global face index of that face
    #   The first NBF entries of global_face_dict are all of the boundary faces
    global_face_dict = {}
    
    fg = 0         # Counts the global faces
    face_hits = {} # Tracks if a face has shown up twice. If it hasn't, this must be a boundary face
    
    # Get the nodes in a given element, organized by element type
    for element_type in gmsh.model.mesh.getElementTypes(3):
        _, nodes = gmsh.model.mesh.getElementsByType(element_type)
        NK = gmsh.model.mesh.getElementProperties(element_type)[3]
        
        # Now define the faces
        for i in range(0, len(nodes), NK):
            element_nodes = nodes[i:i+NK]
            sorted_faces  = get_faces (element_nodes, element_type, sort=True)
            
            # Add faces to the dictionary if they're not already there
            for face in sorted_faces:
                if face not in global_face_dict:
                    fg += 1
                    global_face_dict[face] = fg
                
                if face not in face_hits:
                    face_hits[face] = 0
                
                face_hits[face] += 1
    
    # Now reorganize so that boundary faces are first
    boundary_faces_set = set([face for face, count in face_hits.items() if count == 1])
    sorted_faces_dict = {}
    
    # Add the boundary faces
    fg_boundary = 0
    for face in boundary_faces_set:
        fg_boundary += 1
        sorted_faces_dict[face] = fg_boundary
    
    # Add the interior faces
    fg = fg_boundary
    for face, old_index in global_face_dict.items():
        if face not in boundary_faces_set:
            fg += 1
            sorted_faces_dict[face] = fg
    
    global_face_dict = sorted_faces_dict
    NBF = fg_boundary
    
    return global_face_dict, NBF

def all_faces (global_face_dict, NBF):
    import sys
    # ef2fg is a list of dictionaries. The eth dictionary gives the local face indices with
    # global face indices
    # ELEMENT INDEXING STARTS AT 0
    ef2fg  = []
    
    # fg2els is a dictionary giving the one or two elements sharing a global face
    fg2els = {}
    
    # NKf is an array containing the number of nodes in a global face
    NKf = {}
    
    # Initialize shift to 0, this will be iteratively added to so that
    #   element types other than the first don't get shifted all the way back to zero
    SHIFT = 0
    # Get the nodes in a given element, organized by element type
    for element_type in gmsh.model.mesh.getElementTypes(3):
        element_tags, nodes = gmsh.model.mesh.getElementsByType(element_type)
        NK = gmsh.model.mesh.getElementProperties(element_type)[3]
        
        # Shift the element tags
        SHIFT = SHIFT + min(element_tags) - 1
        
        element_tags = [e - SHIFT for e in element_tags]
        
        # Now define the faces
        for i in range(0, len(nodes), NK):
            element_nodes = nodes[i:i+NK]
            element       = element_tags[i // NK]
            sorted_faces  = get_faces (element_nodes, element_type, sort=True)
            
            lcl_dict = {}
            #print(element_nodes)
            for lcl_face_idx, face_nodes in enumerate(sorted_faces):
                #print(lcl_face_idx + 1, face_nodes)
                fg = global_face_dict[face_nodes]
                #print(fg)
                lcl_dict[lcl_face_idx + 1] = fg
                
                if fg not in NKf:
                    NKf[fg] = len(face_nodes)
                
                if fg not in fg2els:
                    fg2els[fg] = []
                
                fg2els[fg].append(element)
            
            ef2fg.append(lcl_dict)
    
    # This gives the nodes from fg
    # This is basically the inverse of global_face_dict
    fg2kg = {}
    for fg, face_nodes in enumerate(global_face_dict):
        fg2kg[fg + 1] = face_nodes
    
    # This gives the local face index from the element and global face index
    efg2f = {}
    for e, lcl_dict in enumerate(ef2fg):
        efg2f[e + 1] = {fg: f for f, fg in lcl_dict.items()}
    #print(efg2f[100])
    # Now determine the neighbors
    neighbors = {}
    
    for fg, elements in fg2els.items():
        if len(elements) == 1:
            # If there is only one element, this is a boundary element
            e1 = elements[0]
            
            # If this element hasn't been added to neighbors, add it
            if e1 not in neighbors:
                neighbors[e1] = {}
            
            # Now determine the local face index and assign 0 to this neighbor spot
            f = efg2f[e1].get(fg, None) # MAP FROM fg AND e1 TO flocal
            neighbors[e1][f] = 0 # Boundary face (no neighbor)
        
        elif len(elements) == 2:
            # If there are two elements, they must be made to be neighbors
            e1, e2 = elements
            #print(fg, e1, e2)
            # If these elements haven't been added to neighbors, add them
            if e1 not in neighbors:
                neighbors[e1] = {}
            if e2 not in neighbors:
                neighbors[e2] = {}
            
            # Determine local face indices to make assignments
            f = efg2f[e1].get(fg, None) # MAP FROM fg AND e1 TO flocal
            neighbors[e1][f] = e2
            
            f = efg2f[e2].get(fg, None) # MAP FROM fg AND e2 TO flocal
            neighbors[e2][f] = e1
    
    # LOCAL FACES
    fname = os.path.join(bpath, f'{scratch}TAPE%0_local_faces.txt')
    with open(fname, 'w') as f:
        for e, lcl_dict in enumerate(ef2fg):
            fcstr = ' '.join(map(str, list(lcl_dict.values())))
            f.write(f'{e + 1} {fcstr}\n')
    
    # GLOBAL FACES
    fname = os.path.join(bpath, f'{scratch}TAPE%0_global_faces.txt')
    with open(fname, 'w') as f:
        # FORMAT:
        # First entry is the number of faces and number of boundary faces
        # Afterwards, in every odd numbered row, the first column is the 
        # number of nodes in the face, and the second and third columns are
        # the two elements sharing that face (or 0 if its a boundary face, i.e.,
        # the first NBF global faces by my convention)
        # In every even numbered row, it's the global nodes in the face,
        # sorted by size - not according to a particular convention
        f.write(f'{len(fg2els)} {NBF}\n')
        for fg in sorted(fg2els.keys()):
            NK       = NKf[fg]
            elements = fg2els[fg]
            nodes    = fg2kg[fg]
            if len(elements) == 1:
                fcstr    = ' '.join(map(str, elements))
                nodesstr = ' '.join(map(str, nodes))
                f.write(f'{NK} {fcstr} 0\n')
                f.write(f'{nodesstr}\n')
            elif len(elements) == 2:
                fcstr = ' '.join(map(str, elements))
                nodesstr = ' '.join(map(str, nodes))
                f.write(f'{NK} {fcstr}\n')
                f.write(f'{nodesstr}\n')
    
    # NEIGHBORS
    print(f'{scratch}TAPE%0_neighbors.txt')
    fname = os.path.join(bpath, f'{scratch}TAPE%0_neighbors.txt')
    with open(fname, 'w') as f:
        for e in sorted(neighbors.keys()):
            neighstr = ''
            
            for fl, ep in sorted(neighbors[e].items()):
                neighstr = neighstr + ' ' + str(ep)
            
            f.write(f'{e}{neighstr}\n')

debug, bpath, fname, scratch = parser ()

gmsh.initialize()
if (debug):
    gmsh.option.setNumber('General.Terminal', 1)
else:
    gmsh.option.setNumber('General.Terminal', 0)

gmsh.open(fname)

raw_nodes ()

elements ()

global_face_dict, NBF = define_global_faces ()

all_faces (global_face_dict, NBF)

gmsh.finalize()