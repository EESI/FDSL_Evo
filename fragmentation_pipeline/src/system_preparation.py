import os

from biopandas.pdb import PandasPdb

from src.conversion import Converter


class PDBQTprep:
    def __init__(self, protein_file: str, flex_file: str, plip: bool):
        self.prot = protein_file
        self.flex = flex_file
        self.protname = os.path.basename(self.prot)
        self.flexname = os.path.basename(self.flex)
        self.plip = plip

    def files_preparation(self):
        def prepare_ligand():
            flex_df = PandasPdb().read_pdb(self.flex)
            flex_df.df['ATOM'].drop_duplicates(subset='atom_number', keep='first', inplace=True)  # extract first model
            flex_df.df['ATOM']['segment_id'].replace(r'.{1,}', '', regex=True, inplace=True)
            flex_df.df['ATOM']['blank_4'].replace(r'.{1,}', '', regex=True, inplace=True)  # clean pdbqt inheritance
            flex_df.df['ATOM'] = flex_df.df['ATOM'][flex_df.df['ATOM']['element_symbol'] != 'H']  # deprotonate
            return flex_df.to_pdb(path=os.getcwd() + '/{}-prep.pdb'.format(self.flexname.split('.')[0]),
                                  records=['ATOM', 'HETATM'],
                                  gz=False,
                                  append_newline=True)

        def prepare_protein():
            prot_df = PandasPdb().read_pdb(self.prot)
            prot_df.df['ATOM'] = prot_df.df['ATOM'][prot_df.df['ATOM']['element_symbol'] != 'H']  # deprotonate

            prot_df.df['ATOM']['residue_name'].replace(to_replace = 'CYX', value = 'CYS', inplace = True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='HIE', value='HIS', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='HID', value='HIS', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='HIP', value='HIS', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='ASX', value='ASN', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='GLX', value='GLN', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='GLH', value='GLN', inplace=True)
            prot_df.df['ATOM']['residue_name'].replace(to_replace='LYN', value='LYS', inplace=True)

            return prot_df.to_pdb(path=os.getcwd() + '/{}-dehyd.pdb'.format(self.protname.split('.')[0]),
                                  records=['ATOM', 'HETATM'],
                                  gz=False,
                                  append_newline=True)

        prepare_protein()
        prepare_ligand()
        print(self.plip)
        return Converter(os.getcwd() + '/{}-dehyd.pdb'.format(self.protname.split('.')[0]),
                         os.getcwd() + '/{}-prep.pdb'.format(self.flexname.split('.')[0]), self.plip).convert()
