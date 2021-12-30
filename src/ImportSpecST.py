# encoding: iso-8859-15
import bns,shelve
import os
import os.path
import re
import shutil

#Le dialog de base duquel dérive les autres class
Dlg = bns.Windows.XmlDlgBuilder
MessageBox = bns.Util.Program.MessageBox

#Repertoire de base
#BaseDir=r"U:\TEST_BIONUMERICS\SCRIPT\Script_test\\" # le repertoire avec les lots d'analyse SBT
BaseDir = r"S:\partage\BioNumericsAnalyse\LegioSBT\SequencesBrutes\\"
sbt_gene = ['flaA','pilE','asd','mip','mompS','proA','neuA'] # liste des genes SBT

class MainDlg(Dlg.Dialogs):
    """
    Class permettant d'importer dans BioNumerics les resultats SBT enregitres dans le fichier ST_LotID.txt obtenu a partir du site
    http://www.hpa-bioinformatics.org.uk/legionella/legionella_sbt/php/sbt_homepage.php
    """

    def __init__(self):
        Dlg.Dialogs.__init__(self,'MainDlg',remembersettings='local')

        self.list_exp=["flaA","proA","asd","mip","mompS","neuA-neuAh","pilE"]

        #Liste des lots avec fichiers multifasta sbt prets pour importation
        self.listDir=[x for x in os.listdir(BaseDir) if (x.startswith('_') and not x.endswith('-'))]
        #Lot choisis pas l'utilisateur
        self.targetDir=""
        #Lot apparaissant par defaut dans le drop list
        try:
            self.defaultdir=self.listDir[0]
        except:
            self.defaultdir='NA'

        #Sous titre apparaissant dans le dialog
        self.subtitle = Dlg.StaticText(u"Importation des résultats SBT dans BioNumerics",maxWidth=50)
        #Message a l'utilisateur
        self.message = Dlg.StaticText("Veuillez sélectionnez un lot d'analyse: ")

        #La drop list contenant la liste des lots
        self.ListDirDrop = Dlg.Drop("dropList",self.listDir,10,30,canEdit=False,default=self.defaultdir)

        #Le bouton permettant d'amorcer l'importation des resultats SBT
        self.ImportBtn = Dlg.ButtonAction(u"Démarrer l'importation",self.LaunchImport,buttonID="startimport")

        #Disposition des elements dans le dialog
        self.grid = [[self.subtitle],[self.message],[self.ListDirDrop],[self.ImportBtn]]

        #Une instance de la class de base Dlg.SimpleDialog
        self.simpleDlg = Dlg.SimpleDialog(self.grid)

        #Creer le dialog
        self.AddSimpleDialog(u"Legionella SBT - Importation de résultats",self.simpleDlg)

    def FindAlleleSeq(self,exp,allele):
        """
        Retourne la sequence nucleotidique d un allele pour un gene donne
        """

        seq="" # la sequence nucleotidique
        line = "temp" #  contenu temporaire de ligne

        # pour differencier les numeros d allele neuAh de neuA
        if exp=='neuA-neuAh':
            if int(allele) > 200:
                exp='neuAh'
            else:
                exp='neuA'

        with open(r"S:\partage\BioNumericsAnalyse\LegioSBT\ST_ref\ST_seq_corrected.fasta") as readf: # fichier les sequences nucleotidiques des differents alleles

            line = "temp"
            while not (re.search(r'_'+exp+r'_'+allele+'$',line)) and line != '' : # on lie le fichier tant que l on a pas trouve une sequence qui match
                line = readf.readline()

            if re.search(r'_'+exp+r'_'+allele+'$',line):
                #print "line is ", line
                seq=readf.readline() # la sequence nucleotidique

            else:
                print 'No seq allele find for ', exp, " allele ",allele
                seq=""

        readf.close()

        return seq # on retourne la sequence

    def LockSeqByExp(self,key,exp):
        """
        Lock nucleotidic sequence for specific entry key - experience
        """

        if not bns.ConnectedDb.ObjectID(bns.ConnectedDb.ObjectType('SEQUENCES'),{'KEY':key,'EXPERIMENT':exp}).IsLocked():
            bns.ConnectedDb.ObjectID(bns.ConnectedDb.ObjectType('SEQUENCES'),{'KEY':key,'EXPERIMENT':exp}).SetLocked(1)

    def CopyToCons(self):
        """
        Copier les alleles de references dans les alleles correspondantes de specimens cliniques.
        Car les sequences consensus des specimens cliniques de l ancienne bd n ont pas ete trimmer.
        """

        allele=""
        expfield=""

        Entries=bns.Database.Db.Entries
        for entry in Entries:
            if entry.Selection ==True: # les entrees selectionnees
                #print entry.Key
                for exp in self.list_exp:

                    if re.search('LEGIOSBT',bns.Database.Db.Info.Name.upper()) and exp=='neuA-neuAh': # car dans la bd LegioSBT le champ neuA-neuAh se nomme neuA dans le schema MySQL
                        expfield="neuA"
                    else:
                        expfield=exp

                    allele=bns.Database.EntryField(entry.Key, expfield).Content # recuperer l allele du gene pour ce specimen

                    if allele !="": # s il y a une allele

                        myseq=bns.Sequences.Sequence()
                        d = {'seq': ""} # sequence nucleotidique
                        allele_seq=self.FindAlleleSeq(exp,allele) # recuperer la sequence correspondante pour cette allele du gene

                        try: # on cee l exp_cons si elle n existe pas deja
                            myseq.Create(exp+'_cons',"",entry.Key)
                        except:
                            myseq.Load(entry.Key,exp+"_cons") # si l experience existe deja, on la charge
                            #print 'Exp already exist for ', entry.Key

                        myseq.Set(allele_seq) # on set la sequence nucleotidique recupere du fichier ST_seq_corrected.fasta dans l exp_cons
                        myseq.Save()

                        self.LockSeqByExp(entry.Key,exp+'_cons') # on lock la sequence

                        del(myseq)

                    else:
                        print "No allele for ", entry.Key, " in exp ", exp


    def CopyNeuAToNeuANeuAh(self):
        '''
        Copier les sequences de neuA ou neuAh dans neuA-neuAh. OBSOLETE
        '''
        Entries=bns.Database.Db.Entries
        for entry in Entries:
            if entry.Selection ==True: # les entrees selectionnees

                d = {'seq': ""} # sequence nucleotidique

                myseq=bns.Sequences.Sequence()

                if myseq.Load(entry.Key,'neuAh'): # tenter de loader l experience neuAh pour cette key
                    pass

                else: # si neuAh n existe pas pour cette key, on load neuA
                    myseq.Load(entry.Key,'neuA')

                if myseq.Get(byref=d): # recuperer la sequence nucleotidique de l experience

                    try: # on copie la sequence dans l experience neuA-neuAh
                        myseq.Create('neuA-neuAh',"",entry.Key)
                        myseq.Set(d["seq"])
                        myseq.Save()

                    except:
                        print 'Problem with key ', entry.Key

    def LaunchImport(self,args):

        bns.Database.Db.Selection.Clear() # deselection de tous les entrees de la bd

        self.targetDir=self.ListDirDrop.GetValue() # le lot choisis par l utilisateur

        with open(BaseDir+self.targetDir+r"\ST"+self.targetDir+".txt") as readf: # le fichier resultat ST_LotID.txt
           for line in readf:
               try:
                   lineformat = re.search(r'^(\S+)\t(\d+|NA)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)$',line) # format de ligne avec resultat SBT

                   bns.Database.Entry(lineformat.group(1)).Selection = True # selection des specimens du lot

                   bns.Database.EntryField(lineformat.group(1), "Sequence Type").Content = lineformat.group(2) # inscrire le resultat ST dans le champ Sequence Type

                   group = 3
                   for field in sbt_gene: # inscrire les resultats d allele dans les champs appropries

                       if group == 9: # pour inscrire le numero d allele dans le champ neuA ou neuAh. Si > 200 => neuAh

                           bns.Database.EntryField(lineformat.group(1), 'neuA').Content = lineformat.group(group)

                           ''' Ce code n est plus necessaire. On place aussi les neuAh dans le champ neuA
                           if int(lineformat.group(group)) > 200: # allele neuAh
                               bns.Database.EntryField(lineformat.group(1), 'neuAh').Content = lineformat.group(group)
                           else: # allele neuA
                               bns.Database.EntryField(lineformat.group(1), 'neuA').Content = lineformat.group(group)
                           '''

                       else:
                           bns.Database.EntryField(lineformat.group(1), field).Content = lineformat.group(group)
                           group += 1
               except:
                   pass

        readf.close()

        bns.Database.Db.Fields.Save()

        shutil.move(BaseDir+self.targetDir+r"/",BaseDir+self.targetDir+"-/") # ajouter un - devant le nom du repertoire du lot

        #self.CopyNeuAToNeuANeuAh() # copie des sequences consensus des specimens du lot dans l experience neuA-neuAh
        self.CopyToCons() # copie des sequences allele de reference dans les exp_cons correspondantes

        MessageBox(" " ,u'Importation terminée','information')

#Afficher le MainDlg
dlg = MainDlg()
dlg.Show()
