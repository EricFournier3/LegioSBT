# encoding: iso-8859-15
#//$menu window=main;popup=Analysis;name=Template-SBT;
import bns,shelve
import os
import os.path
import re

"""
CDSLSPQ-2120 Eric Fournier 2019-08-12
    - ajout de l'amorce flaA-R-N
    - augmenter de 12 a 13 la longueur maximale de noms de specimen dans les fichiers ab1
"""
 
#Le dialog de base duquel dérive les autres class
Dlg = bns.Windows.XmlDlgBuilder
MessageBox = bns.Util.Program.MessageBox

#Repertoire de base
#rawSeqBaseDir=r"U:\TEST_BIONUMERICS\LegioSBT\\"
rawSeqBaseDir=r"S:\partage\BioNumericsAnalyse\LegioSBT\SequencesBrutes\\"


class MainDlg(Dlg.Dialogs):	
    def __init__(self):

        Dlg.Dialogs.__init__(self,'MainDlg',remembersettings='local')

        #Dictionnary qui associe un extension d'amorce à une expérience Bionumerics pour IDBM
        #CDSLSPQ-2120
        self.expCode={'asd511F':'asd', 'asd1039R':'asd','flaA587F':'flaA','flaA960R':'flaA','mip74F':'mip','mip595R':'mip'
        ,'mompS450F':'mompS','mompS1015R':'mompS','mompS509F':'mompS','neuA196F':'neuA','neuA611R':'neuA','neuAh_L':'neuAh', 'neuAh_R':'neuAh','pilE35F':'pilE','pilE453R':'pilE','proA1107F':'proA','proA1553R':'proA','flaA-R-N':'flaA'}

        #Liste des repertoires avec des sequences ABI bruts
        self.listDir=[]
        #Repertoire choisis pas l'utilisateur
        self.targetDir=""
        #Repertoire apparaissant par defaut dans le drop list
        self.defautdir=""
        self.defaultLot=''

        ListDir = [x for x in os.listdir(rawSeqBaseDir) if (os.path.isdir(rawSeqBaseDir+"\\"+x) and not (str(x).startswith('_')))]
        print ListDir
        ListDirCut=[int(x[0:8]) for x in ListDir] # Supprimer ce qui suit les -
        LotDict=dict(zip(ListDir,ListDirCut))  # associe le nom de lot coupe a son nom original => cree un dict
        LotTupleList=LotDict.items()  # convertir en liste de pair de tuple
        LotTupleSorted=sorted(LotTupleList,key=lambda lot:lot[1],reverse=True) # faire le trie du lot le plus recent au plus ancien

        LotTupleSortedCut=LotTupleSorted[0:6] # conserver les 6 lots les plus recents

        # Ajouter les 6 lots au menu deroulant
        for item in reversed(LotTupleSortedCut):
            self.listDir.append(item[0])
            self.defaultLot = item[0] # le plus recent lot affiche

        #Sous titre apparaissant dans le dialog
        self.subtitle = Dlg.StaticText("Création du template d'assemblage à partir des traces AB1")
        #Message a l'utilisateur
        self.message = Dlg.StaticText("Veuillez sélectionnez le répertoire des traces: ")
        #La drop list contenant la liste des repertoires
        self.ListDirDrop = Dlg.Drop("dropListReperoire",self.listDir,10,30,canEdit=False,default=self.defaultLot)
        #Le bouton permettant de creer le fichier template d'assemblage
        self.CreateTemplateBtn = Dlg.ButtonAction("Créer le template d'assemblage",self.CreateTemplateFile,buttonID="CreateTemplateBouton")
        #Disposition des elements dans le dialog
        self.grid = [[self.subtitle],[self.message],[self.ListDirDrop],[self.CreateTemplateBtn]]

        #Une instance de la class de base Dlg.SimpleDialog
        self.simpleDlg = Dlg.SimpleDialog(self.grid)
        #Creer le dialog self.simpleDlg ayant pour titre IDBM - Template d'assemblage
        self.AddSimpleDialog("Legionella SBT - Template d'assemblage",self.simpleDlg)

    #Creation du fichier template d'assemblage
    def CreateTemplateFile(self,args):
        #Liste avec 2 elements: Le numero du specimen (8 caracteres) et le nom de l'amorce
        IdPrimer=[]
        self.targetDir=self.ListDirDrop.GetValue()

        with open(rawSeqBaseDir+'Template_'+self.targetDir+'.txt', 'w') as f:
            #Le header dans le fichier template
            f.write('[ID]\t[KEY]\t[EXP]\n')
            #Pour tous les fichiers dans le repertoire choisis
            for item in os.listdir(rawSeqBaseDir+self.targetDir):
                 #Si est un fichier .ab1
                 if item.endswith('.ab1'):
                    item_corrected=re.sub('__','_',item) # car certain fichier contiennent un _ supplementaire. Car erreur lors de la preparation de la plaque excel
                    IdPrimer=self.ParseTraceName(item_corrected)
                    #CDSLSPQ-2120
                    parsedItem=re.search(r'(^.{2,13})_(neuAh_[LR]|.*)_\w{3}_\d{4}-\d{2}-\d{2}_\d{1}',item_corrected)
                    #f.write(parsedItem.group()+'\t'+IdPrimer[0]+'\t'+self.expCode[IdPrimer[1]]+'\t'+'\n')
                    f.write(re.sub('\.ab1','',item)+'\t'+IdPrimer[0]+'\t'+self.expCode[IdPrimer[1]]+'\t'+'\n')
        f.close()

        with open(rawSeqBaseDir+self.targetDir+'\Fasta_'+self.targetDir+'.fasta', 'w') as f: # le fichier fasta vide pour contenir les contigs apres edition
            pass
        f.close()

        MessageBox(" " ,'Le template a été créé avec succès','information')

    #Permet d'associer le nom du specimen a une amorce
    def ParseTraceName(self,trace):
        IdPrimer=[]
        #CDSLSPQ-2120
        traceformat=re.search(r'(^.{2,13})_(neuAh_[LR]|.*)_\w{3}_\d{4}-\d{2}-\d{2}',trace)
        if traceformat:
                IdPrimer.append(traceformat.group(1))
                IdPrimer.append(traceformat.group(2))
                #print IdPrimer

        return IdPrimer

    def TestDirPrint(self,args):
        print "targetdir is ", self.targetDir

#Afficher le MainDlg		
dlg = MainDlg()
dlg.Show()

