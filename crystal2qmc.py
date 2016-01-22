import numpy as np
import sys

periodic_table = [
  "h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar",
  "k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br",
  "kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te",
  "i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm",
  "yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
  "fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
  "rf","db","sg","bh","hs","mt","ds","rg","cp","uut","uuq","uup","uuh","uus","uuo"
]

# Reads in the geometry, basis, and pseudopotential from GRED.DAT.
# TODO pseudopotential may not yet work with more than one atom type.
def read_gred():
  lat_parm = {}
  ions = {}
  basis = {}
  pseudo = {}

  gred = open("GRED.DAT",'r').read()

  # Fix numbers with no space between them.
  gred = gred.replace("-"," -")
  gred = gred.replace("E -","E-") 

  gred_words = gred.split()
  nparms = [int(w) for w in gred_words[1:4]]
  cursor = 4

  # These follow naming of cryapi_inp ("inf" -> "info").
  info = [int(w) for w in gred_words[cursor          :cursor+nparms[0]]]
  itol = [int(w) for w in gred_words[cursor+nparms[0]:cursor+nparms[1]]]
  par  = [int(w) for w in gred_words[cursor+nparms[1]:cursor+nparms[2]]]
  cursor += sum(nparms)

  # Lattice parameters.
  lat_parm['prim_cell'] = \
      np.array([float(w) for w in gred_words[cursor:cursor+9]]).reshape(3,3)
  cursor += 9
  prim_trans= np.array([float(w) for w in gred_words[cursor:cursor+9]]).reshape(3,3)
  cursor += 9
  lat_parm['conv_cell'] = prim_trans.dot(lat_parm['prim_cell'])
  cursor += info[1] + 48*48 + 9*info[1] + 3*info[1] # Skip symmetry part.

  # Lattice "stars" (?) skipped.
  cursor += info[4]+1 + info[78]*3 + info[4]+1 + info[4]+1 + info[78] + info[78]*3

  # Some of ion information.
  natoms = info[23]
  ions['charges'] = [float(w) for w in gred_words[cursor:cursor+natoms]]
  cursor += natoms
  # Atom positions.
  atom_poss = np.array([float(w) for w in gred_words[cursor:cursor+3*natoms]])
  ions['positions'] = atom_poss.reshape(natoms,3)
  cursor += 3*natoms

  # Basis information (some ion information mixed in).
  nshells = info[19]
  nprim   = info[74]
  # Formal charge of shell.
  basis['charges'] = np.array([float(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells
  # "Adjoined gaussian" of shells.
  basis['adj_gaus'] = np.array([float(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells
  # Position of shell.
  shell_poss = np.array([float(w) for w in gred_words[cursor:cursor+3*nshells]])
  basis['positions'] = shell_poss.reshape(nshells,3)
  cursor += 3*nshells
  # Primitive gaussian exponents.
  basis['prim_gaus'] = np.array([float(w) for w in gred_words[cursor:cursor+nprim]])
  cursor += nprim
  # Coefficients of s, p, d, and (?).
  basis['coef_s'] = np.array([float(w) for w in gred_words[cursor:cursor+nprim]])
  cursor += nprim
  basis['coef_p'] = np.array([float(w) for w in gred_words[cursor:cursor+nprim]])
  cursor += nprim
  basis['coef_dfg'] = np.array([float(w) for w in gred_words[cursor:cursor+nprim]])
  cursor += nprim
  basis['coef_max'] = np.array([float(w) for w in gred_words[cursor:cursor+nprim]])
  cursor += nprim
  # Skip "old normalization"
  cursor += 2*nprim
  # Atomic numbers.
  ions['atom_nums'] = np.array([int(w) for w in gred_words[cursor:cursor+natoms]])
  cursor += natoms
  # First shell of each atom (skip extra number after).
  basis['first_shell'] = np.array([int(w) for w in gred_words[cursor:cursor+natoms]])
  cursor += natoms + 1
  # First primitive of each shell (skips an extra number after).
  basis['first_prim'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells + 1
  # Number of prims per shell.
  basis['prim_shell'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells
  # Type of shell: 0=s,1=sp,2=p,3=d,4=f.
  basis['shell_type'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells
  # Number of atomic orbtials per shell.
  basis['nao_shell'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells
  # First atomic orbtial per shell (skip extra number after).
  basis['first_ao'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells + 1
  # Atom to which each shell belongs.
  basis['atom_shell'] = np.array([int(w) for w in gred_words[cursor:cursor+nshells]])
  cursor += nshells

  # Pseudopotential information.
  # Skipping what I assume is what type of pseudopot input it is.
  cursor += natoms+1
  ngauss = int(gred_words[cursor])
  cursor += 1
  headlen = int(gred_words[cursor])
  cursor += 1
  # (?)
  itpse = int(gred_words[cursor])
  cursor += 1
  # Exponents of r^l prefactor.
  pseudo['r_exps'] = -1*np.array([int(w) for w in gred_words[cursor:cursor+ngauss]])
  cursor += ngauss
  # Number of Gaussians for angular momenutum j
  pseudo['n_per_j'] = np.array([int(w) for w in gred_words[cursor:cursor+headlen]])
  cursor += headlen
  # (?) Not used.
  nsom = [int(w) for w in gred_words[cursor:cursor+itpse]]
  cursor += itpse + 1
  # Actual floats of pseudopotential.
  pseudo['exponents'] = np.array([float(w) for w in gred_words[cursor:cursor+ngauss]])
  cursor += ngauss
  pseudo['prefactors'] = np.array([float(w) for w in gred_words[cursor:cursor+ngauss]])
  cursor += ngauss

  return info, lat_parm, ions, basis, pseudo

# Reads in kpoints and eigenvectors from KRED.DAT.
def read_kred(info,basis):
  eigsys = {}

  kred = open("KRED.DAT",'r').read()
  kred_words = kred.split()
  cursor = 0

  # Number of k-points in each direction.
  eigsys['nkpts_dir'] = [int(w) for w in kred_words[cursor:cursor+3]]
  cursor += 3
  # Total number of inequivilent k-points.
  nikpts = int(kred_words[cursor])
  cursor += 1
  # Reciprocal basis.
  recip_vecs = np.array([float(w) for w in kred_words[cursor:cursor+9]])
  eigsys['recip_vecs'] = recip_vecs.reshape(3,3)
  cursor += 9
  # Inequivilent k-point coord in reciprocal basis.
  ikpt_coords = np.array([float(w) for w in kred_words[cursor:cursor+3*nikpts]])
  ikpt_coords = ikpt_coords.reshape(nikpts,3)
  cursor += 3*nikpts
  # is complex (0) or not (1), converted to True (if complex) or False
  eigsys['ikpt_iscmpx'] = \
    np.array([int(w) for w in kred_words[cursor:cursor+nikpts]]) == 0
  cursor += nikpts
  # Skip symmetry information.
  cursor += 9*48
  # Geometric weight of kpoints.
  eigsys['kpt_weights'] = [float(w) for w in kred_words[cursor:cursor+nikpts]]
  cursor += nikpts
  # Eigenvalues: (how many) = (spin) * (number of basis) * (number of kpoints)
  nevals = (info[63]+1)*info[6]*nikpts
  eigsys['eigvals'] = [float(w) for w in kred_words[cursor:cursor+nevals]]
  cursor += nevals
  # Weights of eigenvales--incorperating Fermi energy cutoff.
  eig_weights = [float(w) for w in kred_words[cursor:cursor+nevals]]
  cursor += nevals

  # Read in eigenvectors at inequivilent kpoints. Can't do all kpoints because we 
  # don't know if non-inequivilent kpoints are real or complex (without symmetry
  # info)
  nevals_kpt = int(round(nevals / nikpts))
  nkpts  = np.prod(eigsys['nkpts_dir'])
  nao = sum(basis['nao_shell'])
  ncpnts = nevals_kpt * nao
  nreads = 0
  kpt_coords = []
  eigvecs    = []
  for kidx in range(nkpts):
    #kpt_coords[kidx,:] = np.array([int(w) for w in kred_words[cursor:cursor+3]])
    new_kpt_coord = np.array([int(w) for w in kred_words[cursor:cursor+3]])
    cursor += 3

    # If new_kpt_coord is an inequivilent point...
    ikidx = np.prod(new_kpt_coord == ikpt_coords,1) == 1
    if any(ikidx):
      if eigsys['ikpt_iscmpx'][ikidx]: nreads = 2*ncpnts
      else:                            nreads = ncpnts
      eig_k = np.array([float(w) for w in kred_words[cursor:cursor+nreads]])
      cursor += nreads
      kpt_coords.append(new_kpt_coord)
      eigvecs.append(eig_k.reshape(int(round(nreads/nao)),nao))
    else: # ...else, skip.
      skip = True
      while skip:
        try: # If there's an int, we're at next kpoint.
          int(kred_words[cursor])
          skip = False
        except ValueError: # Keep skipping.
          cursor += ncpnts
        except IndexError: # End of file.
          skip = False
          break

  # It's probably true that eigsys['kpt_coords'] == ikpt_coords,
  # because we only read in inequivilent kpoints. However, ordering might be
  # different, and the ordering is correct for eigsys['kpt_coords'].
  eigsys['kpt_coords'] = kpt_coords
  eigsys['eigvecs'] = eigvecs

  return eigsys

def find_basis_cutoff(lat_parm):
  latvec = lat_parm['prim_cell']
  cutoff_divider = 2.000001
  cross01 = np.cross(latvec[0], latvec[1])
  cross12 = np.cross(latvec[1], latvec[2])
  cross02 = np.cross(latvec[0], latvec[2])

  heights = [0,0,0]
  heights[0]=abs(np.dot(latvec[0], cross12)/np.dot(cross12,cross12)**.5)
  heights[1]=abs(np.dot(latvec[1], cross02)/np.dot(cross02,cross02)**.5)
  heights[2]=abs(np.dot(latvec[2], cross01)/np.dot(cross01,cross01)**.5)
  return min(heights)/cutoff_divider

# TODO generalize to ferro.
# TODO generalize to spin-polarized.
def write_slater(basis,kidx,base="qwalk"):
  kbase = base + '_' + str(kidx)
  nmo = sum(basis['charges']) / 2
  if nmo % 1 > 1e-10: print("Error: number of electrons is probably noninteger")
  nmo = int(round(nmo))
  outlines = [
      "SLATER",
      "ORBITALS {",
      "CUTOFF_MO",
      "  MAGNIFY 1",
      "  NMO {0}".format(nmo),
      "  ORBFILE {0}.orb".format(kbase),
      "  INCLUDE {0}.basis".format(base),
      "  CENTERS { USEGLOBAL }",
      "}",
      "DETWT { 1.0 }",
      "STATES {",
      "  # Spin up orbitals.", 
      "  " + " ".join([str(i+1) for i in range(nmo)]),
      "  # Spin down orbitals.",
      "  " + " ".join([str(i+1) for i in range(nmo)]),
      "}"
    ]
  with open(kbase+".slater",'w') as outf:
    outf.write("\n".join(outlines))
  return outlines # Might be confusing.

# f orbital normalizations are from 
# <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
def normalize_eigvec(eigsys,basis,kidx):
  eigvec = eigsys['eigvecs'][kidx]
  snorm = 1./(4.*np.pi)**0.5
  pnorm = snorm*(3.)**.5
  dnorms = [
      .5*(5./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      .5*(15./(4.*np.pi))**.5,
      (15./(4*np.pi))**.5
    ]
  fnorms = [
      ( 7./(16.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (105./(16.*np.pi))**.5, # xyz
      (105./(4.*np.pi))**.5,
      (35./(32.*np.pi))**.5,
      (35./(32.*np.pi))**.5
    ]

  # Duplicate coefficients for complex, and if multiple basis elements are d.
  # This is to align properly with the d-components of eigvecs.
  if eigsys['ikpt_iscmpx'][kidx]:
    tmp = [[f,f] for f in dnorms]
    dnorms = []
    for l in tmp: dnorms += l
  tmp = [[f for i in range(sum(basis['shell_type']==3))] for f in dnorms]
  dnorms = []
  for l in tmp: dnorms += l
  dnorms = np.array(dnorms)
  # Likewise for f.
  if eigsys['ikpt_iscmpx'][kidx]:
    tmp = [[f,f] for f in fnorms]
    fnorms = []
    for l in tmp: fnorms += l
  tmp = [[f for i in range(sum(basis['shell_type']==4))] for f in fnorms]
  fnorms = []
  for l in tmp: fnorms += l
  fnorms = np.array(fnorms)

  ao_type = []
  for sidx in range(len(basis['shell_type'])):
    if eigsys['ikpt_iscmpx'][kidx]:
      ao_type += \
        [basis['shell_type'][sidx] for ao in range(2*basis['nao_shell'][sidx])]
    else:
      ao_type += \
        [basis['shell_type'][sidx] for ao in range(basis['nao_shell'][sidx])]
  ao_type = np.array(ao_type)

  if any(ao_type==1):
    print("sp orbtials not implemented in normalize_eigvec(...)")
    exit("Not implemented.")


  eigvec[ao_type==0,:] *= snorm
  eigvec[ao_type==2,:] *= pnorm
  eigvec[ao_type==3,:] = ((eigvec[ao_type==3,:].T)*dnorms).T
  eigvec[ao_type==4,:] = ((eigvec[ao_type==4,:].T)*fnorms).T
  return None
      
def write_orb(eigsys,basis,ions,kidx,base="qwalk"):
  outf = open(base + '_' + str(kidx) + ".orb",'w')
  eigvecs = eigsys['eigvecs'][kidx]
  nao_atom = int(sum(basis['nao_shell']) / len(ions['positions']))
  coef_cnt = 1
  for moidx in np.arange(eigvecs.shape[1])+1:
    for atidx in np.unique(basis['atom_shell']):
      for aoidx in np.arange(nao_atom)+1:
        outf.write(" ".join(map(str,[moidx,aoidx,atidx,coef_cnt]))+"\n")
        coef_cnt += 1
  coef_cnt -= 1 # Last increment doesn't count.
  if eigsys['ikpt_iscmpx'][kidx]: ncoefs = eigvecs.size//2
  else:                           ncoefs = eigvecs.size
  if coef_cnt != ncoefs:
    print("Error: Number of coefficients not coming out correctly!")
    print("Counted: {0} \nAvailable: {1}".format(coef_cnt,eigvecs.size))
    sys.exit("Debug Error")
  eigvecs_flat = eigvecs.flatten()
  print_cnt = 0
  outf.write("COEFFICIENTS\n")
  for cidx in range(coef_cnt):
    if eigsys['ikpt_iscmpx'][kidx]:
      outf.write("({0},{1}) ".format(eigvecs_flat[2*cidx],eigvecs_flat[2*cidx+1]))
    else:
      outf.write("{0} ".format(eigvecs_flat[cidx]))
    print_cnt += 1
    if print_cnt % 5 == 0: outf.write("\n")
  outf.close()
  return None

# TODO Generalize to ferro.
# TODO Generalize to spin-polarized.
# TODO Molecule.
# TODO Generalize to no pseudopotential.
# TODO pseudopotential may not yet work with more than one atom type.
# TODO Don't understand ordering of pseudopotential.
def write_sys(lat_parm,basis,eigsys,pseudo,kidx,base="qwalk"):
  min_exp = min(basis['prim_gaus'])
  cutoff_length = (-np.log(1e-8)/min_exp)**.5
  basis_cutoff = find_basis_cutoff(lat_parm)
  cutoff_divider = basis_cutoff*2.0 / cutoff_length
  kbase = base + '_' + str(kidx)
  nmo = sum(basis['charges']) / 2
  if nmo % 1 > 1e-10: print("Error: number of electrons is probably noninteger")
  nmo = int(round(nmo))
  outlines = [
      "system { periodic",
      "  nspin {{ {0} {1} }}".format(nmo,nmo),
      "  latticevec {",
    ]
  for i in range(3):
    outlines.append("    " + " ".join(map(str,lat_parm['prim_cell'][i])))
  outlines += [
      "  }",
      "  origin { 0 0 0 }",
      "  cutoff_divider {0}".format(cutoff_divider),
      "  kpoint {{ {0} }}".format(" ".join(map(str,eigsys['kpt_coords'][kidx])))
    ]
  for aidx in range(len(ions['positions'])):
    outlines.append(
      "  atom {{ {0} {1} coor {2} }}".format(
        periodic_table[ions['atom_nums'][aidx]-200-1], # Assumes ECP.
        ions['charges'][aidx],
        " ".join(map(str,ions['positions'][aidx]))
      )
    )
  outlines.append("}")
  for aidx in np.unique(ions['atom_nums']):
    atom_name = periodic_table[aidx-200-1]
    atom_name = periodic_table[aidx-200-1]
    numL = sum(pseudo['n_per_j']>0)

    for i in range(1,len(pseudo['n_per_j'])):
      if (pseudo['n_per_j'][i-1]==0)and(pseudo['n_per_j'][i]!=0):
        print("Weird pseudopotential, please generalize write_sys(...).")
        sys.exit("Not implemented.")

    n_per_j = pseudo['n_per_j'][pseudo['n_per_j']>0]
    order = list(range(n_per_j[0],sum(n_per_j)))+list(range(n_per_j[0]))
    exponents   = pseudo['exponents'][order]
    prefactors  = pseudo['prefactors'][order]
    r_exps      = pseudo['r_exps'][order]
    if numL > 2: aip = 12
    else:        aip =  6
    outlines += [
        "pseudo {",
        "  {0}".format(atom_name),
        "  aip {0}".format(aip),
        "  basis {{ {0}".format(atom_name),
        "    rgaussian",
        "    oldqmc {",
        "      0.0 {0}".format(numL),
        "      "+" ".join(map(str,n_per_j[1:]))+" "+str(n_per_j[0])
      ]
    cnt = 0
    for jidx,j in enumerate(n_per_j):
      for eidx in range(j):
        outlines.append("      " + " ".join(map(str,[
          r_exps[cnt]+2,
          exponents[cnt],
          prefactors[cnt]
        ])))
        cnt += 1
    outlines += ["    }","  }","}"]
  with open(kbase+".sys",'w') as outf:
    outf.write("\n".join(outlines))
  return None

# TODO check generality of numbers.
# TODO generalize atom_type.
def write_jast2(lat_parm,base="qwalk"):
  basis_cutoff = find_basis_cutoff(lat_parm)
  atom_type = "si"
  outlines = [
      "jastrow2",
      "group {",
      "  optimizebasis",
      "  eebasis {",
      "    ee",
      "    cutoff_cusp",
      "    gamma 24.0",
      "    cusp 1.0",
      "    cutoff {0}".format(basis_cutoff),
      "  }",
      "  eebasis {",
      "    ee",
      "    cutoff_cusp",
      "    gamma 24.0",
      "    cusp 1.0",
      "    cutoff {0}".format(basis_cutoff),
      "  }",
      "  twobody_spin {",
      "    freeze",
      "    like_coefficients { 0.25 0.0 }",
      "    unlike_coefficients { 0.0 0.5 }",
      "  }",
      "}",
      "group {",
      "  optimize_basis",
      "  eibasis {",
      "    {0}".format(atom_type),
      "    polypade",
      "    beta0 0.2",
      "    nfunc 3",
      "    rcut {0}".format(basis_cutoff),
      "  }",
      "  onebody {",
      "    coefficients {{ {0} 0.0 0.0 0.0}}".format(atom_type),
      "  }",
      "  eebasis {",
      "    ee",
      "    polypade",
      "    beta0 0.5",
      "    nfunc 3",
      "    rcut {0}".format(basis_cutoff),
      "  }",
      "  twobody {",
      "    coefficients { 0.0 0.0 0.0 }",
      "  }",
      "}"
    ]
  with open(base+".jast2",'w') as outf:
    outf.write("\n".join(outlines))
  return None

# TODO Test for more than one atom type.
def write_basis(basis,ions,base="qwalk"):
  hybridized_check = 0.0
  hybridized_check += sum(abs(basis['coef_s'] * basis['coef_p']))
  hybridized_check += sum(abs(basis['coef_p'] * basis['coef_dfg']))
  hybridized_check += sum(abs(basis['coef_s'] * basis['coef_dfg']))
  if hybridized_check > 1e-10:
    print("Hybridized AOs (like sp) not implmemented in write_basis(...)")
    sys.exit("Not implemented.")

  # If there's no hybridization, at most one of coef_s, coef_p, and coef_dfg is
  # nonzero. Just add them, so we have one array.
  done_atoms = {}
  coefs = basis['coef_s'] + basis['coef_p'] + basis['coef_dfg']

  shell_type = np.tile("Unknown...",basis['shell_type'].shape)
  typemap = ["S","SP","P","5D","7F_crystal","G","H"]
  for i in range(5): shell_type[basis['shell_type']==i] = typemap[i]

  cnt = 0
  aidx = 0
  atom_type = ions['atom_nums'][aidx]
  outlines = [
      "basis {",
      "  {0}".format(periodic_table[atom_type-200-1]),
      "  aospline",
      "  normtype CRYSTAL",
      "  gamess {"
    ]
  for sidx in range(len(shell_type)):
    new_aidx = basis['atom_shell'][sidx]-1

    new_atom_type = ions['atom_nums'][new_aidx]
    if aidx != new_aidx:
      done_atoms[atom_type] = True
      if new_atom_type in done_atoms.keys(): 
        continue
      else:
        atom_type = new_atom_type
        aidx = new_aidx
        outlines += [
            "basis {",
            "  {0}".format(periodic_table[atom_type-200-1]),
            "  aospline",
            "  normtype CRYSTAL",
            "  gamess {"
          ]

    nprim = basis['prim_shell'][sidx]
    outlines.append("    {0} {1}".format(shell_type[sidx],nprim))
    for pidx in range(nprim):
      outlines.append("      {0} {1} {2}".format(
        pidx+1,
        basis['prim_gaus'][cnt],
        coefs[cnt]
      ))
      cnt += 1
  outlines += ["  }","}"]
  with open(base+".basis",'w') as outf:
    outf.write("\n".join(outlines))
  return None

def write_moanalysis():
  return None

#########################
# Begin actual execution.
info, lat_parm, ions, basis, pseudo = read_gred()
eigsys = read_kred(info,basis)
for (kidx,kpt) in enumerate(eigsys['kpt_coords']):
  write_slater(basis,kidx)
  normalize_eigvec(eigsys,basis,kidx)
  write_orb(eigsys,basis,ions,kidx)
  write_sys(lat_parm,basis,eigsys,pseudo,kidx)
  write_basis(basis,ions)
  write_jast2(lat_parm)
