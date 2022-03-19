import os
import freud
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib

def getGroups(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail):
    """Get the groups involved in Hydrogen bonding.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atoms

    Returns
    -------

    OxygenHead : List of integres
        List of acceptor oxygen atoms.
    Hydrogen : List of integers
        List of donor Hydrogen atoms.
    OxygenTail : List of integers
        List of oxygen atoms bonded to the donor H atom.
    list_names_hydrogen : List of strings
        List of the names of the donor H atoms.
    list_names_oxygen_head : List of strings
        List of the names of the acceptor O atoms.
    list_names_oxygen_tail : List of strings
        List of the names of the O atoms boned to donor H atoms.
    """



    top = traj.topology
    OxygenHead = top.select(sel_oxygen_head)
    Hydrogen = top.select(sel_hydrogen)
    OxygenTail = top.select(sel_oxygen_tail)

    return OxygenHead, Hydrogen, OxygenTail, list_names_hydrogen,list_names_oxygen_head, list_names_oxygen_tail 

def ellipticalFun(dist,cosangle):
    """Check if a point falls in the elliptical region.
    Parameters
    ----------
    dist : float
        Distance in nm
    cosangle : float
        cosine of the angle

    Returns
    -------
    Boolean
        Returns True if the point lies in the elliptical region. Otherwise False.
    """

    if ((((dist-0.30)/0.05)**2 + (1/0.3572**2)*(cosangle+1)**2) < 1): return True

    else: return False

def create_bond_dict(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided=False, max_r_bond=0.12):
    """If the bond information is not provided in the traj.top file, then this function helps create bonds between hydrogens and oxygens based on a distance cutoff.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atomsi
    bonded_pdb_provided : Boolean
        Boolean containing the info if a PDB with bond information is provided on not
    max_r_bond : float
        Distance (nm) cutoff for bond between O and H

    Returns
    -------
    OxygenTail_bond_diction : Dictionary
        Dictory containing the information about bond between atoms so that the information is required to be retrieved from the top inforamation in the main loop
    top : mdTraj topology
        Updated topology after bond information updation
    """




    top = traj.top
    atoms = list(top.atoms)
    OxygenHead, Hydrogen, OxygenTail, Hydrogen_atom_names, OxygenHead_names, OxygenTail_names = getGroups(traj,sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
    OxygenTail_bond_diction={}
    
    if bonded_pdb_provided == False:
        for OxygenTail_index in  OxygenTail:
            
            OxygenTail_bond_diction[OxygenTail_index]=[]
            
        frame = 0
        L = traj[frame].unitcell_lengths
        #print(L)
        box = freud.box.Box(Lx=L[0][0],Ly= L[0][1], Lz=L[0][2])
        points_O_tail=traj.xyz[frame][OxygenTail]
        points_H = traj.xyz[frame][Hydrogen]
        aq = freud.locality.AABBQuery(box, points_H)


        query_result = aq.query(points_O_tail, dict(r_max=max_r_bond))
        for bond in query_result:
            if np.isclose(bond[2], 0): #means that the distance is zero
                continue
            point = bond[0]
            neighbor = bond[1]
            dist = bond[2]
            #print(point, neighbor)
            original_Otail_index = OxygenTail[point]
            original_H_index = Hydrogen[neighbor]
            OxygenTail_bond_diction[original_Otail_index].append(original_H_index)
            top.add_bond(atoms[original_Otail_index],atoms[original_H_index ])

                    
                    
        return OxygenTail_bond_diction, top

        

    for OxygenTail_index in  OxygenTail:
        print("first_loop")
        OxygenTail_bond_diction[OxygenTail_index]=[]
        for bond in top.bonds:
            print(bond.atom1.name)
            print(bond.atom2.name)
            if (bond.atom1.index == OxygenTail_index) or (bond.atom2.index == OxygenTail_index):

                if bond.atom1.name in Hydrogen_atom_names:
                    bonded_H_index = bond.atom1.index

                elif bond.atom2.name in Hydrogen_atom_names:
                    bonded_H_index = bond.atom2.index

                else:
                    continue
                OxygenTail_bond_diction[OxygenTail_index].append(bonded_H_index)

            else:
                continue
    return OxygenTail_bond_diction


def calualateHBMap(traj, r_cutoff, res_r, res_a, skip_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail):
    """Function for calculating the 2D histrogram (cos angle vs distance) for determing H bonds.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    r_cutoff : float
        Distance (nm) up to which pairs are considered.
    res_r : int
        Number of bins in the r direction
    res_a : int
        Number of bins in the cos angle direction
    skip_frames : int
        Number of frames to be skipped at the start of the traj (maybe for eqlb).
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atomsi

    Returns
    -------
    rdf_output : numpy array
        bins for r (nm)
    inter_output : numpy array
        bins for cos angle
    map_output : numpy 2D array
        2D array for the frequency on the 2D grid
    ans : ndarray
        n_bonds x 2 array that contains the pair of indices of atoms that make a H bond
    hbond_time : array
        number of Hbonds in each frame


    Examples
    --------
    >>> res_r = 200 ; res_a = 200 ; r_cutoff = 0.75 ; skip =0
    >>> sel_oxygen_head = 'name O5' ; sel_hydrogen = 'name H1 or name H2' ; sel_oxygen_tail = 'name O5'; list_names_hydrogen = ["H1", "H2"] ; list_names_oxygen_head = ["O5"] ; list_names_oxygen_tail = ["O5"]
    >>> rdf_output, inter_output, map_output,hbond,hbond_time = calualateHBMap(traj, r_cutoff, res_r, res_a, skip, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
    >>> #Plotting
    >>> plt.figure() ; cmap = plt.get_cmap('jet') ; plt.figure(figsize=(5, 3)) ; plt.style.use('default'); levels = np.linspace(0,10,11) ; cs = plt.contourf(rdf_output[0], inter_output[0], map_output,levels = levels, cmap=cmap) ; plt.xlabel('r (nm)') ; plt.ylabel('cos(\u03B8)') ; plt.xlim([0.2, 0.4]) ; plt.ylim([-1, -0.5]) ; plt.colorbar()

    """

    bond_diction, top = create_bond_dict(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
    traj.top = top
    top = traj.top
    print("traj as {} bonds".format(top.n_bonds))
    dr = r_cutoff/res_r
#    radii = np.linspace(0.0, res_r * dr, res_r)

    radii = [dr/2]
    for i in range(1,res_r):
        radii.append(radii[i-1]+dr)
    radii = np.array(radii)


    angles = np.linspace(-1, 1, num=res_a)

    da = angles[1]-angles[0]
    print("da", da)
    print("dr", dr)
    volume = np.zeros(res_r);
    for i in range(0, len(radii)-1):
        volume[i] = (4/3)*np.pi*(radii[i+1]**3-radii[i]**3)
    volume[-1]=((4/3)*np.pi*((radii[-1]+dr)**3-radii[-1]**3))

    """ Initialize output array """
    gr_final = np.zeros(res_r)
    intangle_final = np.zeros(res_a)
    maps_final = np.zeros((res_a,res_r))
    hbond_time = []
    ans = []
    prev = 0

    """ Iterate through each frame """
    for frame in range(skip_frames, traj.n_frames):

        print("working on {}".format(frame))

        L = traj[frame].unitcell_lengths
        #print(L)
        box = freud.box.Box(Lx=L[0][0],Ly= L[0][1], Lz=L[0][2])
        intangle_prob = np.zeros(res_a)
        g_of_r = np.zeros(res_r)
        maps = np.zeros((res_a,res_r))

        """ Get Neighbors """
        OxygenHead, Hydrogen, OxygenTail, Hydrogen_atom_names, OxygenHead_names, OxygenTail_names = getGroups(traj,  sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
        numHead = len(OxygenHead)
        numTail = len(OxygenTail)
        points_O_tail=traj.xyz[frame][OxygenTail]
        points_O_head=traj.xyz[frame][OxygenHead]
        aq = freud.locality.AABBQuery(box, points_O_head)



        """ Iterate through Neighbors """

        query_result = aq.query(points_O_tail, dict(r_max=r_cutoff))
        pairs_analyzed = []

        for shbond in query_result:
            if np.isclose(shbond[2], 0): #means that the distance is zero
                continue
            point = shbond[0]
            neighbor = shbond[1]
            dist = shbond[2]
            #print(point, neighbor)
            original_Otail_index = OxygenTail[point]
            original_Ohead_index = OxygenHead[neighbor]
            if [original_Otail_index, original_Ohead_index] in pairs_analyzed:
                continue
            else:
                pairs_analyzed.append([original_Otail_index, original_Ohead_index])
                pairs_analyzed.append([original_Ohead_index, original_Otail_index])
            #now check if tail is connected to a H
            #print(bond_diction)
            for bonded_H in bond_diction[original_Otail_index]:

                original_H_index = bonded_H
               # print("{}-{} is bonded to {}-{}".format(original_Otail_index,bond.atom1.name, original_H_index,bond.atom2.name))
                index_d = int(dist / dr)
                if 0 < index_d < res_r:
                    g_of_r[index_d] += 1.00
                    """ calculate angle """
                    normal_angle = md.compute_angles(traj[frame], np.array([[original_Otail_index, original_H_index, original_Ohead_index]]))[0][0]

                    inter_angle = np.cos(normal_angle) #np.cos(normal_angle/180*np.pi)
                    #print(inter_angle)
                    if inter_angle < 1 and inter_angle > -1:
                        index_i = np.digitize(inter_angle,angles)-1
                        intangle_prob[index_i] += 1.00
                        """ put into map """
                        maps[index_i,index_d] += 1.00
                        if ellipticalFun(dist, inter_angle) :
                            ans.append([original_Otail_index, original_Ohead_index])


        hbond_time.append(len(ans)-prev)
        prev = len(ans)
        """ normalization """
        total_g = numHead * (numTail-1) / (L[0][0]*L[0][1]*L[0][2])
        #total_g = numHead * (numTail) / (L[0][0]*L[0][1]*L[0][2])
        total_m = np.zeros(res_a)

        """ map normalization """
        for i, value in enumerate(g_of_r):
            g_of_r[i] = value / volume[i] /total_g

        for i in range(0, res_a):

            for j in range(0, res_r):

                maps[i,j] /=  volume[j] * total_g *da



        #g_of_r /= sum(g_of_r)*dr
        intangle_prob /= da*sum(intangle_prob)

        intangle_final += intangle_prob
        maps_final += maps
        gr_final += g_of_r

    intangle_final /= frame
    maps_final /= frame
    gr_final /= frame

    inter_output = np.array([angles, intangle_final])
    rdf_output = np.array([radii, gr_final])
    map_output = maps_final
    return rdf_output, inter_output, map_output,ans,hbond_time


