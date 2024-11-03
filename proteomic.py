#####
# Dr Hugo Boisaubert
# hugo.boisaubert@univ-nantes.fr
# D'après Germain Salvato-Vallverdu (germain.vallverdu@univ-pau.fr) : https://gvallverdu.gitbooks.io/python_sciences/content/genetique.html
#####

acideAmine = {
    "Alanine": {
        "A": "A",
        "Abr": "Ala",
        "masse": 89.09404,
        "pI": 6.00,
        "polaire": False
    },
    "Arginine": {
        "A": "R",
        "Abr": "Arg",
        "masse": 174.20274,
        "pI": 10.76,
        "polaire": True
    },
    "Asparagine": {
        "A": "N",
        "Abr": "Asn",
        "masse": 132.11904,
        "pI": 5.41,
        "polaire": True
    },
    "Aspartate": {
        "A": "D",
        "Abr": "Asp",
        "masse": 133.10384,
        "pI": 2.77,
        "polaire": True
    },
    "Cystéine": {
        "A": "C",
        "Abr": "Cys",
        "masse": 121.15404,
        "pI": 5.07,
        "polaire": False
    },
    "Glutamate": {
        "A": "E",
        "Abr": "Glu",
        "masse": 147.13074,
        "pI": 3.22,
        "polaire": True
    },
    "Glutamine": {
        "A": "Q",
        "Abr": "Gln",
        "masse": 146.14594,
        "pI": 5.65,
        "polaire": True
    },
    "Glycine": {
        "A": "G",
        "Abr": "Gly",
        "masse": 75.06714,
        "pI": 5.97,
        "polaire": False
    },
    "Histidine": {
        "A": "H",
        "Abr": "His",
        "masse": 155.15634,
        "pI": 7.59,
        "polaire": True
    },
    "Isoleucine": {
        "A": "I",
        "Abr": "Ile",
        "masse": 131.17464,
        "pI": 6.02,
        "polaire": False
    },
    "Leucine": {
        "A": "L",
        "Abr": "Leu",
        "masse": 131.17464,
        "pI": 5.98,
        "polaire": False
    },
    "Lysine": {
        "A": "K",
        "Abr": "Lys",
        "masse": 146.18934,
        "pI": 9.74,
        "polaire": True
    },
    "Methionine": {
        "A": "M",
        "Abr": "Met",
        "masse": 149.20784,
        "pI": 5.74,
        "polaire": False
    },
    "Phénylalanine": {
        "A": "F",
        "Abr": "Phe",
        "masse": 165.19184,
        "pI": 5.48,
        "polaire": False
    },
    "Proline": {
        "A": "P",
        "Abr": "Pro",
        "masse": 115.13194,
        "pI": 6.30,
        "polaire": False
    },
    "Sérine": {
        "A": "S",
        "Abr": "Ser",
        "masse": 105.09344,
        "pI": 5.68,
        "polaire": True
    },
    "Thréonine": {
        "A": "T",
        "Abr": "Thr",
        "masse": 119.12034,
        "pI": 5.60,
        "polaire": True
    },
    "Tryptophane": {
        "A": "W",
        "Abr": "Trp",
        "masse": 204.22844,
        "pI": 5.89,
        "polaire": False
    },
    "Tyrosine": {
        "A": "Y",
        "Abr": "Tyr",
        "masse": 181.19124,
        "pI": 5.66,
        "polaire": True
    },
    "Valine": {
        "A": "V",
        "Abr": "Val",
        "masse": 117.14784,
        "pI": 5.96,
        "polaire": False
    }
}


geneticCode = {"Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
               "Phe": ["UUU", "UUC"],
               "Ile": ["AUU", "AUC", "AUA"],
               "Met": ["AUG"],
               "Val": ["GUU", "GUC", "GUA", "GUG"],
               "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
               "Pro": ["CCU", "CCC", "CCA", "CCG"],
               "Thr": ["ACU", "ACC", "ACA", "ACG"],
               "Ala": ["GCU", "GCC", "GCA", "GCG"],
               "Tyr": ["UAU", "UAC"],
               "STOP": ["UAA", "UAG", "UGA"],
               "His": ["CAU", "CAC"],
               "Gln": ["CAA", "CAG"],
               "Asn": ["AAU", "AAC"],
               "Lys": ["AAA", "AAG"],
               "Asp": ["GAU", "GAC"],
               "Glu": ["GAA", "GAG"],
               "Cys": ["UGU", "UGC"],
               "Trp": ["UGG"],
               "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
               "Gly": ["GGU", "GGC", "GGA", "GGG"]}

def codon_to_aa(uncodon:str)-> str:
    """ Renvoie le code à trois lettres d'un acide aminé correspondant au codon
    donné en argument.

    Args:
        uncodon (str): Codon dont on veux l'acide aminé correspondant

    Raises:
        ValueError: Si le codon est inconnu

    Returns:
        str: code à trois lettres de l'acide aminé correspondant
    """

      
    acideAmine = None
    identify = False
    for aa, codons in geneticCode.items():
        if uncodon in codons:
            acideAmine = aa
            identify = True
            break
            
    if not identify:
        raise ValueError("ERREUR : codon '%s' non identifié" % uncodon)

    return acideAmine

def traduction(fragment:str)->str:
    """Traduit le brin d'ARN en séquence peptidique.

    Args:
        fragment (str): séquence à traduire

    Returns:
        str: séquence peptidique produite
    """

    #nombre de codons dans le fragment
    ncodon = len(fragment) // 3
    
    # traduction    
    sequence = ""
    n = 0
    while n < ncodon:
        aa = codon_to_aa(fragment[3*n : 3*n+3])
        if aa != "STOP":
            sequence += aa + "-"
            n += 1
        else:
            sequence += aa
            break
    return sequence

def get_nombre_aa(sequence:str)->int:
    """Retourne le nombre d'acides aminés dans la séquence. Le codon STOP n'est
    pas compté comme un acide aminé. 

    Args:
        sequence (str): sequence à analyser

    Returns:
        int: nombre d'acides aminés dans la séquence
    """
    return len(sequence.split("-")[:-1])

def get_polarite(sequence:str)->float:
    """Retourne le pourcentage d'acides aminés polaires dans une séquence.

    Args:
        sequence (str): sequence à analyser

    Returns:
        float: pourcentage d'acides aminés polaires
    """
    # liste des acides aminés sauf codon STOP
    listeaa = sequence.split("-")[:-1]
    
    npolaire = 0
    for aa in listeaa:
        for aaname in acideAmine:
            if acideAmine[aaname]["Abr"] == aa and acideAmine[aaname]["polaire"]:
                npolaire += 1
    return npolaire / len(listeaa) * 100.

def get_stat_aa(sequence:str)->dict:
    """Compte le nombre de chaque type d'acide aminé et retourne un dictionnaire

    Args:
        sequence (str): sequence à analyser

    Returns:
        dict: nombre de chaque type d'acide aminé
    """
    # liste des acides aminés sauf codon STOP
    listeaa_sequence = sequence.split("-")[:-1]
    
    # liste des 20 codes a trois lettres des acides aminés
    listeaa = [acideAmine[aa]["Abr"] for aa in acideAmine]
    
    # statistique
    data = dict()
    for aa in listeaa:
        naa = listeaa_sequence.count(aa)
        if naa != 0:
            data[aa] = naa
            data[aa] = 0
    
    return data
