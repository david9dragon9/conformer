from .aligngap import align
from .py_qcprot.py_qcprot import rmsd

tto = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def extract_ca(pdb_file, target_chain_id):
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    results = []
    for line in lines:
        if line[:4] == "ATOM":
            comp, num, atype, rtype, chainid, rnum, x, y, z = (
                line[:4],
                line[6:11],
                line[12:16],
                line[16:20],
                line[21],
                line[22:26],
                line[30:38],
                line[38:46],
                line[46:54],
            )
            comp, num, atype, rtype, chainid, rnum, x, y, z = (
                comp.strip(),
                num.strip(),
                atype.strip(),
                rtype.strip(),
                chainid.strip(),
                rnum.strip(),
                x.strip(),
                y.strip(),
                z.strip(),
            )
            num, rnum, x, y, z = int(num), int(rnum), float(x), float(y), float(z)
            if atype == "CA" and chainid == target_chain_id and rtype in tto:
                results.append([rtype, x, y, z])
    return results


def smart_rmsd(pdb1, pdb2, chain1="A", chain2="A", print_stuff=True):
    protein1 = extract_ca(pdb1, chain1)
    protein2 = extract_ca(pdb2, chain2)
    if print_stuff:
        print(pdb1, pdb2, chain1, chain2)
        print(f"Protein 1 {len(protein1)}, Protein 2 {len(protein2)}", flush=True)
    if len(protein1) < 2000 and len(protein2) < 2000:
        one1 = "".join([tto[protein1[i][0]] for i in range(len(protein1))])
        one2 = "".join([tto[protein2[i][0]] for i in range(len(protein2))])
        aligned1, aligned2 = align(one1, one2)
        if print_stuff:
            print(f"Aligned Residues: {len(aligned1)}")
        sc1 = [protein1[x][1:] for x in aligned1]
        sc2 = [protein2[x][1:] for x in aligned2]
        return rmsd(sc1, sc2)
    else:
        if print_stuff:
            print("Skipped due to length.")
