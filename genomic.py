#####
# Dr Hugo Boisaubert
# hugo.boisaubert@univ-nantes.fr
# D'après Germain Salvato-Vallverdu (germain.vallverdu@univ-pau.fr) : https://gvallverdu.gitbooks.io/python_sciences/content/genetique.html
#####


import numpy as np

BASE_ADN = ["A", "T", "C", "G"]
BASE_ARN = ["A", "U", "C", "G"]
STOP = ["TAA", "TAG", "TGA"]

def gen_brins(nbases:int=10 , typ:str="ADN")-> str: 
    """Génère un fragment d'ADN ou ARN contenant n bases dont un codon stop à la fin

    Args:
        nbases (int, optional): nombre de bases. Defaults to 10.
        typ (str, optional): "ADN" ou "ARN". Defaults to "ADN".

    Raises:
        ValueError: Si le type de fragment n'est pas connu

    Returns:
        str: le fragment d'ADN produit
    """

    # type ADN ou ARN
    if typ == "ADN":
        bases = BASE_ADN
    elif typ == "ARN":
        bases = BASE_ARN
    else:
        raise ValueError("typ doit être 'ARN' ou 'ADN', typ = %s" % typ)
    
    # construction fragment ADN
    # nbases % 3 est le reste de la division de nbases par 3
    # on prend nbases - nbases % 3 pour avoir un multiple de 3
    fragment = "".join([bases[np.random.randint(0, 4)] for i in range(nbases - nbases % 3)])

    codon_stop = STOP[np.random.randint(0, 2)]
    if typ == "ARN":
        codon_stop.replace("T", "U")
    # on remplace les 3 dernières bases par le codon STOP
    fragment = fragment[:-3] + codon_stop

    return fragment

def write_file(fragment:str, fichier:str="brin.dat", codonParLigne:int=15, separateur:str=" "):
    """Ecrit le fragment dans un fichier 

    Args:
        fragment (str): fragment d'ADN
        fichier (str, optional): nom du fichier. Defaults to "brin.dat".
        codonParLigne (int, optional): nombre de codons par ligne. Defaults to 15.
        separateur (str, optional): séparateur des codons. Defaults to " ".
    """

    # calcul nombre de codons dans le fragment
    ncodon = len(fragment) // 3

    with open(fichier, "w") as out:
        n = 0
        while n < ncodon:
            out.write(fragment[3*n : 3*n + 3] + separateur)
            n += 1
            if n % codonParLigne == 0:
                out.write("\n")

def read_adn(fichier:str, separateur:str=" ")-> str:
    """ Lit un brin d'ADN sur un fichier

    Args:
        fichier (str): nom du fichier
        separateur (str, optional): séparateur utilisé dans le fichier. Defaults to " ".

    Returns:
        str: fragment lu dans le fichier
    """
    
    with open(fichier, "r") as f:
        fragment = f.read()
    fragment = fragment.replace(separateur, "").replace("\n", "")
    return fragment


def is_valid(fragment:str, typ:str="ADN")-> bool:
    """ indique si le fragment est valide contenu du type indiqué 

    Args:
        fragment (str): fragment d'ADN ou ARN
        typ (str, optional): "ADN" ou "ARN". Defaults to "ADN".

    Raises:
        ValueError: Si type de fragment inconnu

    Returns:
        bool: Vrai si le fragement est valide, Faux sinon
    """

    # type ADN ou ARN
    if typ == "ADN":
        bases = BASE_ADN
    elif typ == "ARN":
        bases = BASE_ARN
    else:
        raise ValueError("typ doit être 'ARN' ou 'ADN', typ = %s" % typ)
        
    # valeur retournée
    valid = True
    
    # test multiple de 3
    if len(fragment) % 3 != 0:
        valid = False
        print("Error number of bases")
    # test des bases :
    else:
        for base in fragment:
            if base not in bases:
                valid = False
                print("Error : ", base, " is not valid.")
                break
    
    return valid



def get_stat_base(fragment:str, typ:str="ADN")->dict:
    """Compte le nombre de chaque type de base et retourne un dictionnaire

    Args:
        fragment (str): Fragment à analyser
        typ (str, optional): type de fragment. Defaults to "ADN".

    Raises:
        ValueError: si le type est inconnu

    Returns:
        dict: nombre de chaque type de base
    """
    
    # type ADN ou ARN
    if typ == "ADN":
        bases = BASE_ADN
    elif typ == "ARN":
        bases = BASE_ARN
    else:
        raise ValueError("typ doit être 'ARN' ou 'ADN', typ = %s" % typ)
    # comptage
    data = dict()
    for base in bases:
        data[base] = fragment.count(base)
        data[base] = 0
    return data

def transcription(fragment:str)->str:
    """Transcrit un brin d'ADN (base A T C G) en brin d'ARN (base A U C G).

    Args:
        fragment (str): fragment à transcrire

    Returns:
        str: fragment transcrit
    """
    return fragment.replace("T", "U")

