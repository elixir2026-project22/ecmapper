"""Synonymous BRENDA cofactor labels (different ligand IDs or spellings for the same
compound) -> one label. Used by brenda_map_cofactors.py (and, idempotently, by
merge_roche_cofactors.py). Specific compounds are not merged into generic classes
(UQ/UQH2 stay separate from "a quinone"/"a quinol")."""
BRENDA_SYNONYMS = {
    "reduced NADPH-hemoprotein reductase": "[reduced NADPH-hemoprotein reductase]",
    "reduced [NADPH-hemoprotein reductase]": "[reduced NADPH-hemoprotein reductase]",
    "oxidized NADPH-hemoprotein reductase": "[oxidized NADPH-hemoprotein reductase]",
    "oxidized [NADPH-hemoprotein reductase]": "[oxidized NADPH-hemoprotein reductase]",
    "acceptor": "oxidized electron acceptor", "an acceptor": "oxidized electron acceptor",
    "oxidized acceptor": "oxidized electron acceptor", "electron acceptor": "oxidized electron acceptor",
    "an oxidized electron acceptor": "oxidized electron acceptor", "A": "oxidized electron acceptor",
    "reduced acceptor": "reduced electron acceptor", "a reduced acceptor": "reduced electron acceptor",
    "a reduced electron acceptor": "reduced electron acceptor", "AH2": "reduced electron acceptor",
    "reduced ferredoxin [Fe-S] cluster": "Fd_red", "oxidized ferredoxin [Fe-S] cluster": "Fd_ox",
    "thioredoxin": "red. thioredoxin", "thioredoxin disulfide": "ox. thioredoxin",
    "electron-transfer flavoprotein": "electron-transferring flavoprotein",
    "oxidized electron-transferring flavoprotein": "electron-transferring flavoprotein",
    "oxidized electron transfer flavoprotein": "electron-transferring flavoprotein",
    "reduced electron-transfer flavoprotein": "reduced electron-transferring flavoprotein",
    "3'-phosphoadenylylsulfate": "3'-phosphoadenylyl sulfate",
    "reduced adrenal ferredoxin": "reduced adrenodoxin", "oxidized adrenal ferredoxin": "oxidized adrenodoxin",
    "reduced quinone": "a quinol",
    "UDP-Gal": "UDP-alpha-D-galactose", "UDP-glucuronate": "UDP-alpha-D-glucuronate",
    "UDP-N-acetyl-alpha-D-glucosamine": "UDP-alpha-N-acetyl-D-glucosamine",
    "a carboxylate": "carboxylate",
    "4a-hydroxy-tetrahydrobiopterin": "4a-hydroxytetrahydrobiopterin",
    # label on ligand ID 17 (diphosphate) that names DMAPP
    "dimethylallyl-PP": "dimethylallyl diphosphate",
}
