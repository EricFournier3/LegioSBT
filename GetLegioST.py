# encoding: iso-8859-15
#//$menu window=main;popup=Analysis;name=Get ST;
import bns,shelve
import os
import os.path
import re
from Bio import SeqFeature
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import subprocess

"""
CDSLSPQ-2244 
Eric Fournier 2019-10-03 Creation de ce script pour LegioSBT blast local
"""

#TODO database update

Dlg = bns.Windows.XmlDlgBuilder
MessageBox = bns.Util.Program.MessageBox

#Nom de la base de données BioNumerics
DB_NAME = bns.Database.Db.Info.Name.upper()

#Les noms de champs dans BioNumerics
if re.search('DEV', DB_NAME):
    BION_FIELD_NAME = {'ST':'SEQUENCE_TYPE','FLAA':'FLAA','PILE':'PILE','ASD':'ASD','MIP':'MIP','MOMPS':'MOMPS','PROA':'PROA','NEUA':'NEUA'}
else:
    BION_FIELD_NAME = {'ST': 'Sequence Type', 'FLAA': 'flaA', 'PILE': 'pilE', 'ASD': 'asd', 'MIP': 'mip','MOMPS': 'mompS', 'PROA': 'proA', 'NEUA': 'neuA'}

#Prefix du nom de fichier multi-fasta
SeqInPrefix = "Fasta_"

#suffix des noms de fichier
spec_file_suffix = "_spec.fasta" # fichiers de sequences specifique a chaque locus
db_file_suffix = "_db" # fichiers de base de donnee
blast_file_suffix = "_blast.txt" # fichiers de resultats blast

#pour les path des fichiers selon le locus
sbt_file_dict = {"flaA":{"spec_file":"","db_file":"","blast_file":""},"pilE":{"spec_file":"","db_file":"","blast_file":""},"asd":{"spec_file":"","db_file":"","blast_file":""},
            "mip":{"spec_file":"","db_file":"","blast_file":""},"mompS":{"spec_file":"","db_file":"","blast_file":""},"proA":{"spec_file":"","db_file":"","blast_file":""},
            "neuA":{"spec_file":"","db_file":"","blast_file":""}}

#data pour chacun des specimens
spec_dict = {}

#les records pour chacun des locus
rec_dict = {"flaA":[],"pilE":[],"asd":[],"mip":[],"mompS":[],"proA":[],"neuA":[]}

#Les parametres blast
blast_param = {'BlastTool': r"\\swsfi03p\Partage\BioNumericsAnalyse\BlastToolKit\blast-2.6.0+\bin\blastn",'max_target_seqs':1, 'perc_identity': 100,'outfmt':7}

#pour mode debug
debug = False

#path vers database et lot d'analyse
if debug:
    lot_basedir = os.path.join("U:\\","PROJETS","DT1","C2","CDSLSPQ-2244","SequencesBrutes")
    db_basedir = os.path.join("U:\\","PROJETS","DT1","C2","CDSLSPQ-2244","SBT_DB")
else:
    lot_basedir = os.path.join("S:\\", "Partage", "BioNumericsAnalyse", "LegioSBT","SequencesBrutes")
    db_basedir = os.path.join("S:\\", "Partage", "BioNumericsAnalyse", "LegioSBT", "SBT_DB")

#fichier st de reference
st_ref_file = os.path.join(db_basedir,'legionella_sbt_ref.txt')

#associer les st de reference au numeros d'alleles
st_ref_dict = {}


def SetSBT_File_Dict(workpath):
    """
    Setter les path vers le fichiers du sbt_file_dict pour chacun
    des locus
    :param workpath:
    :return:
    """
    for gene in sbt_file_dict:
        sbt_file_dict[gene]["spec_file"] = os.path.join(workpath,gene + spec_file_suffix)
        sbt_file_dict[gene]["db_file"] = os.path.join(db_basedir, gene + db_file_suffix)
        sbt_file_dict[gene]["blast_file"] = os.path.join(workpath, gene + blast_file_suffix)

def Set_Spec_Dict(fastahandle):
    """
    Mettre les sequences des specimens dans le spec_dict

    :param fastahandle:
    :return:
    """
    for myrec in fastahandle:
        spec_id = myrec.id
        dna_seq = str(myrec.seq)
        gene = re.search(r'^(\S+)\|(\S+)', spec_id).group(1)
        spec_name = re.search(r'^(\S+)\|(\S+)', spec_id).group(2)

        if not spec_dict.has_key(spec_name):
            spec_dict[spec_name] = {}

        if gene == 'neuAh':
            gene = 'neuA'
        spec_dict[spec_name][gene + "_seq"] = dna_seq


def Set_Rec_Dict():
    """
    On construit le rec_dict pour chacun des locus

    :return:
    """
    for spec in spec_dict.keys():
        for gene in spec_dict[spec].keys():
            gene_name = gene[:-4]
            rec = SeqRecord(Seq(str(spec_dict[spec][gene]), IUPAC.ambiguous_dna),id=gene_name + '|' + spec,description='')

            rec_dict[gene_name].append(rec)

def WriteRecToFile():
    """
    On met les sequences dans les fichiers spec_file
    :return:
    """
    for gene in rec_dict.keys():
        SeqIO.write(rec_dict[gene],sbt_file_dict[gene]['spec_file'],'fasta')

def Blast_On_SBT_DB():
    """
    On blast les seuquences de chaque locus
    :return:
    """

    for gene in sbt_file_dict.keys():
        seq_file = sbt_file_dict[gene]['spec_file']
        blast_file = sbt_file_dict[gene]['blast_file']
        db_file = sbt_file_dict[gene]['db_file']

        ret = subprocess.call([blast_param['BlastTool'],'-db',db_file,'-query',seq_file,'-outfmt',str(blast_param['outfmt']),'-max_target_seqs',str(blast_param['max_target_seqs']),'-perc_identity',str(blast_param['perc_identity']),'-out',blast_file],shell=True)

def SetSpecAllele():
    """
    A partir des fichiers blast on recupere les numeros d'alleles de chaque locus pour chaque specimen
    :return:
    """
    gene_list = sbt_file_dict.keys()

    for spec in spec_dict.keys():
        for gene in gene_list:
            spec_dict[spec][gene + '_allele'] = "NA"

    for gene in gene_list:
        for line in open(sbt_file_dict[gene]["blast_file"]):
            if not line.startswith('#'):
                res = line.split('\t')
                spec = re.search(r'\S+\|(\S+)', res[0]).group(1)
                aln_length = int(res[3])
                sub_start = int(res[8])
                sub_stop = int(res[9])
                subjt_length = sub_stop - sub_start + 1
                allele = re.search(r'\S+_(\d+)', res[1]).group(1)

                if aln_length == subjt_length:
                    # print 'allele is ' + allele
                    # print 'gene is ' + gene
                    spec_dict[spec][gene + '_allele'] = int(allele)

def BuildSTrefDict():
    """
    On contruite le dictionnaire ST de reference versus alleles
    :return:
    """
    st_ref_file_handle = open(st_ref_file)
    st_ref_file_handle.readline()
    for line in st_ref_file_handle:
        if not line.startswith('0'):
            data = line.split('\t')
            data[-1] = str(data[-1]).strip('\n')
            st_ref_dict[data[0]] = {}
            st_ref_dict[data[0]]['flaA'] = int(data[1])
            st_ref_dict[data[0]]['pilE'] = int(data[2])
            st_ref_dict[data[0]]['asd'] = int(data[3])
            st_ref_dict[data[0]]['mip'] = int(data[4])
            st_ref_dict[data[0]]['mompS'] = int(data[5])
            st_ref_dict[data[0]]['proA'] = int(data[6])
            st_ref_dict[data[0]]['neuA'] = int(data[7])

    st_ref_file_handle.close()

def CreateAllelesDict(dict_val):
    """
    Obtenir un dictionnaire contenant seulement les numeros d allele pour chacun
    des locus
    :param dict_val:
    :return:
    """
    alleles_dict = {}

    for key in dict_val.keys():
        if key.endswith('_allele'):
            alleles_dict.update({key[:-7]:dict_val[key]})

    return alleles_dict

def SetSpecST():
    """
    On recherche les ST pour chacun des specimens selon les valeurs d'allele
    :return:
    """

    for spec,dict_val in spec_dict.items():
        found = False
        alleles_dict = CreateAllelesDict(dict_val)

        for ref_st, ref_alleles in st_ref_dict.items():
            if alleles_dict == ref_alleles:

                spec_dict[spec]['st'] = int(ref_st)
                found = True
                break

        if not found:
            spec_dict[spec]['st'] = 'NA'

def ImportInBioNumerics():
    """
    Importation des resulats dans BioNumerics
    :return:
    """
    bns.Database.Db.Selection.Clear()

    for spec, dict_val in spec_dict.items():


        bns.Database.EntryField(spec, BION_FIELD_NAME['ST']).Content = str(spec_dict[spec]['st'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['FLAA']).Content = str(spec_dict[spec]['flaA_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['PILE']).Content = str(spec_dict[spec]['pilE_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['ASD']).Content = str(spec_dict[spec]['asd_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['MIP']).Content = str(spec_dict[spec]['mip_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['MOMPS']).Content = str(spec_dict[spec]['mompS_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['PROA']).Content = str(spec_dict[spec]['proA_allele'])
        bns.Database.EntryField(spec, BION_FIELD_NAME['NEUA']).Content = str(spec_dict[spec]['neuA_allele'])

    bns.Database.Db.Fields.Save()

def RenameLot(lot_path):
    import shutil

    new_lot = '_' + os.path.basename(lot_path) + '\\'
    new_path = os.path.join(os.path.dirname(lot_path),new_lot)

    #os.rename(lot_path, new_path)
    #
    print lot_path + '\\'
    print new_path

    try:
        shutil.move(lot_path + '\\', new_path)
        
    except :
        MessageBox("Erreur", "Impossible de renommer le répertoire", 'exclamation')



def Get_SBT_ST(lot):
    lot_path = os.path.join(lot_basedir, lot)
    work_path = os.path.join(lot_path, "WorkDir")
    if not os.path.isdir(work_path):
        os.mkdir(work_path)

    #fichier fasta exporté du BioNumerics
    fasta_file_name = SeqInPrefix + os.path.basename(lot_path) + ".fasta"
    fasta_handle = SeqIO.parse(os.path.join(lot_path, fasta_file_name), 'fasta')
    __bnscontext__.SetBusy("Traitement en cours...")  # Message d attente
    SetSBT_File_Dict(work_path)
    Set_Spec_Dict(fasta_handle)
    Set_Rec_Dict()
    WriteRecToFile()
    Blast_On_SBT_DB()
    SetSpecAllele()
    BuildSTrefDict()
    SetSpecST()
    ImportInBioNumerics()
    #RenameLot(lot_path) bug
    __bnscontext__.SetBusy("")

    MessageBox(" ", 'Terminé', 'information')


class MainDlg(Dlg.Dialogs):
    """
    Pour l'interface client
    """
    def __init__(self):

        Dlg.Dialogs.__init__(self,'MainDlg',remembersettings='local')

        self.listDir = []
        self.SetListLotDir()

        self.defaultLot = ""

        self.subtitle = Dlg.StaticText("")
        self.message = Dlg.StaticText("Veuillez choisir un lot à traiter: ")
        self.ListDirDrop = Dlg.Drop("dropListReperoire", self.listDir, 10, 30, canEdit=False, default=self.defaultLot)
        self.GetSTBtn = Dlg.ButtonAction("Obtenir les ST", self.Get_SBT_ST,buttonID="getSTBouton")
        self.grid = [[self.subtitle], [self.message], [self.ListDirDrop], [self.GetSTBtn]]
        self.simpleDlg = Dlg.SimpleDialog(self.grid)
        self.AddSimpleDialog("Legionella SBT", self.simpleDlg)

    def SetListLotDir(self):
        ListDir = [x for x in os.listdir(lot_basedir) if (os.path.isdir(os.path.join(lot_basedir,x)) and not (str(x).startswith('_')))]
        ListDirCut = [int(x[0:8]) for x in ListDir]
        DictDir = dict(zip(ListDir, ListDirCut))
        TupleDir = DictDir.items()
        TupleDirSorted = sorted(TupleDir, key=lambda lot: lot[1],reverse=True)
        TupleDirSorted = TupleDirSorted[0:6]

        for item in reversed(TupleDirSorted):
            self.listDir.append(item[0])
            self.defaultLot = item[0]

    def Get_SBT_ST(self, args):
        if not len(self.ListDirDrop.GetValue()) == 0:
            Get_SBT_ST(self.ListDirDrop.GetValue())
        else:
            MessageBox("Erreur", "Choisir un lot d'analyse", 'exclamation')


dlg = MainDlg()
dlg.Show()

