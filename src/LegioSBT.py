# encoding: iso-8859-15
__author__ = 'ericfournier'

import re
from Bio import SeqIO
import mechanize
from selenium import webdriver
import time
from ttk import *
import ttk
from Tkinter import*
import tkMessageBox as box
import os
import shutil

class STfileParser():
    """
    Pour extraire les alleles - ST du fichier SBT_ST.txt obtenu du site www.hpa-bioinformatics.org.uk/legionella.
    Single use.
    """

    def __init__(self):

        self.st_dict = {} # key=ST et value=profile allele

        self.Parse() #
        self.sortST()
        self.WriteST()

    def sortST(self):
        '''
        trie selon les valeurs de ST
        '''
        self.st_tuple = self.st_dict.items()

        self.st_tuple = sorted(self.st_tuple,key=lambda x:int(x[0]))
        #print self.st_tuple

    def Parse(self):
        '''
        Extraire les alleles  et ST correspondant du fichier SBT_ST.txt
        '''

        with open(r'/home/ericfournier/TESTING/PYTHON/Exercices/SBT_ST.txt') as readf:
            for line in readf:

                try:

                    Allele_ST_format = re.search(r'(\d+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)\s+(\d+)',line) # format flaA,pilE,asd,mip,mompS,proA,neuA   ST

                    self.st_dict[Allele_ST_format.group(8)] = [Allele_ST_format.group(1),Allele_ST_format.group(2),Allele_ST_format.group(3),Allele_ST_format.group(4),Allele_ST_format.group(5),Allele_ST_format.group(6),
                                                               Allele_ST_format.group(7)]

                except:

                    pass

        readf.close()


    def WriteST(self):
        '''
        Transferer les allele et ST correspondant dans un fichier formatte
        '''

        with open(r'/home/ericfournier/TESTING/PYTHON/Exercices/SBT_ST_parsed.txt','w') as writef:

            writef.write('ST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA\n\n') # header

            for item in self.st_tuple:
                #print item
                try:

                    writef.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(item[0],item[1][0],item[1][1],item[1][2],item[1][3],item[1][4],item[1][5],item[1][6]))

                except IndexError as e:
                    print 'Error ', e
                    break

        writef.close()


class SBT():

    def __init__(self,lotID):


        self.st_dict = {}  # ST - profil d allele

        self.SBTgene = ['flaA','pilE','asd','mip','mompS','proA','neuA','neuAh'] # liste de genes  legionella SBT

        self.specSet = set()  # liste de specimens clinique

        self.records = None  # fichier fasta des sequences nucleotidiques SBT

        self.baseDir=r"S:\partage\BioNumericsAnalyse\LegioSBT\SequencesBrutes\\" # repertoire avec les lots de sequences brutes
        #self.baseDir=r"C:\Users\ericfournier\Documents\PYTHONtest\SequencesBrutes\\" # repertoire avec les lots de sequences brutes

        self.lotID = lotID # le lot de sequences ab1

        self.spec_seq = {} # key=specimen - valeur=ses sequences nucleotidique sous forme de dictionnaire {flaA:seq,pilE:seq,asd:seq,mip:seq,mompS:seq,proA:seq,neuA:seq}
        self.spec_allele = {} # key=specimen - valeur=son profil allelique sous forme de liste {flaA:allele,pilE:allele,asd:allele,mip:allele,mompS:allele,proA:allele,neuA:allele,ST:st} (obtenu a partir du site www.hpa-bioinformatics.org.uk/legionella)



    def ReadST(self):
        '''
        Mettre les ST - profil d alleles dans un dictionnaire
        '''

        with open(r'S:\partage\BioNumericsAnalyse\LegioSBT\ST_ref\SBT_ST_parsed.txt') as readf:
        #with open(r'C:\Users\ericfournier\Documents\PYTHONtest\ST_ref\SBT_ST_parsed.txt') as readf:
            for line in readf:
                try:
                    line_format = re.search(r'^(\d+)\t(\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+)$',line)
                    self.st_dict[int(line_format.group(1))] = line_format.group(2).split('\t')
                except:
                    pass

        readf.close()

    def ReadSeq(self):
        '''
        Parsing des sequences de specimens du fichiers fasta exporte du BioNumerics
        '''
        self.records = list(SeqIO.parse('{0}{1}/Fasta_{1}.fasta'.format(self.baseDir,self.lotID), "fasta")) # les sequences fasta

        for sequences in self.records:
            header_format = re.search(r'\S+\|(\S+)',sequences.id)
            self.specSet.add(header_format.group(1)) # la liste des specimens.

        for spec in self.specSet: # initialisation d un dictionnaire pour chacun des specimens
            self.spec_seq[spec]={}

        for sequences in self.records:
            header_format = re.search(r'(\S+)\|(\S+)',sequences.id) # un header fasta d une sequence
            self.spec_seq[header_format.group(2)][header_format.group(1)] = sequences.seq.__str__() # attribution de la sequence au gene correspondant dans le dictionnaire du specimen.

    def WriteSpecAlleleProfil(self):
        '''
        Sauvegarde des resultats de profil allelique des specimens cliniques dans ST_LotID.txt
        '''

        STdict = {}

        for spec in self.specSet: # on a de besoin du ST du specimen
            STdict[spec] = self.FindST(spec)

        with open(self.baseDir+self.lotID+'\ST_'+self.lotID+'.txt','w') as writef:
            writef.write('Specimen\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA\n\n')

            for spec in self.specSet:
                #print STdict[spec]
                #print '\t'.join([self.spec_allele[spec][gene] for gene in self.SBTgene])
                #print self.spec_allele[spec], " ", STdict[spec]
                writef.write(spec+'\t'+STdict[spec]+'\t'+'\t'.join([self.spec_allele[spec][gene] for gene in self.SBTgene[0:7]])+'\n')

        writef.close()


    def FindST(self,spec):
        '''
        Recherche du ST correspondant dans le fichier SBT_ST_parsed.txt selon le profil allelique du specimen
        '''

        for seqtype in self.st_dict:
                try:
                    if self.st_dict[seqtype] == [self.spec_allele[spec][gene] for gene in self.SBTgene[0:7]]:
                        #print "Spec ",spec, " ST ", seqtype

                        return  str(seqtype)
                        break
                    else:
                        pass

                except:

                    print 'error'

        return 'NA'

    def SaveAllele(self):
        '''
        Pour extraire en fasta toutes les alleles d'une gene
        Single use
        '''


        find_alleles = []

        rec = list(SeqIO.parse('/home/ericfournier/TESTING/PYTHON/Exercices/ST_seq.fasta', "fasta"))

        with open(r'/home/ericfournier/TESTING/PYTHON/Exercices/neuA_alleles.fasta','w') as writef:

            for sequences in rec:
                header_format = re.search(r'ST_\d+_(\S+)_(\d+)$',sequences.id) # un header fasta d une sequence

                if (header_format.group(1) == 'neuA') and (int(header_format.group(2)) not in find_alleles) and (int(header_format.group(2)) < 59):
                    writef.write('>'+sequences.id+'\n'+sequences.seq.__str__()+'\n')
                    find_alleles.append(int(header_format.group(2)))


        writef.close()

        print sorted(find_alleles)

        missing_alleles = []

        for allele in range(1,59):
            if allele not in find_alleles:
                missing_alleles.append(allele)

        print 'Missing alleles ',missing_alleles


    def SplitST_seq(self):
        """
        Creer des fichiers fasta separes pour chacun des genes du fichiers ST_seq_corrected.fasta pour importation dans BioNumerics
        Single use
        """


        handle = {'flaA':open(self.baseDir+"ST_seq_flaA.fasta", "a"),'pilE':open(self.baseDir+"ST_seq_pilE.fasta", "a"),'asd':open(self.baseDir+"ST_seq_asd.fasta", "a"),
                  'mip':open(self.baseDir+"ST_seq_mip.fasta", "a"),'mompS':open(self.baseDir+"ST_seq_mompS.fasta", "a"),'proA':open(self.baseDir+"ST_seq_proA.fasta", "a"),
                  'neuA':open(self.baseDir+"ST_seq_neuA.fasta", "a"),'neuAh':open(self.baseDir+"ST_seq_neuAh.fasta", "a")}


        rec = list(SeqIO.parse('/home/ericfournier/TESTING/PYTHON/Exercices/ST_seq_corrected.fasta', "fasta"))

        for sequences in rec:
                SeqIO.write(sequences,handle[re.search(r'ST_\d+_(\S+)_\d+$',sequences.id).group(1)], "fasta")



    def neuA_to_neuAh(self):
        """
        Modifier le nom neuA pour neuAh pour les alleles superieures a 200 dans les header fasta du fichier ST_seq.fasta
        Single use
        """

        rec = list(SeqIO.parse('/home/ericfournier/TESTING/PYTHON/Exercices/ST_seq.fasta', "fasta"))

        with open(r'/home/ericfournier/TESTING/PYTHON/Exercices/ST_seq_corrected.fasta','w') as writef:

            for sequences in rec:
                header_format = re.search(r'ST_\d+_(\S+)_(\d+)$',sequences.id) # un header fasta d une sequence

                if (header_format.group(1) == 'neuA') and (int(header_format.group(2)) > 200):
                    sequences.id = re.sub(r'neuA',r'neuAh',sequences.id)
                    writef.write('>'+sequences.id+'\n'+sequences.seq.__str__()+'\n')
                    print sequences.id

                else:
                    writef.write('>'+sequences.id+'\n'+sequences.seq.__str__()+'\n')

        writef.close()

    def CreateFastaST(self):
        '''
        Creation du fichier fasta de ST allele a partir de www.hpa-bioinformatics.org.uk/legionella
        pour importation dans BioNumerics.
        Single use.
        '''

        gene_index = {0:'flaA',1:'pilE',2:'asd',3:'mip',4:'mompS',5:'proA',6:'neuA'}

        with open(r'/home/ericfournier/TESTING/PYTHON/Exercices/ST_seq.fasta','a') as writef:
            # Pour chacun des 2145 ST, enregistre en fasta la sequence de chacun des alleles correspondants
            for st in sorted(self.st_dict.keys())[200:202]:

                #print 'ST is ', st
                index=0

                for allele in self.st_dict[st]:
                    writef.write('>ST_'+str(st)+'_'+gene_index[index]+'_'+str(allele)+'\n') # header fasta
                    dna_seq=''
                    wait=3

                    while dna_seq == '': # tant que l on a pas obtenu une sequence
                        if wait == 15:
                            exit() # on arrete le script si le temps d attente devient trop long
                        driver = webdriver.Firefox()
                        driver.get("http://www.hpa-bioinformatics.org.uk/legionella/legionella_sbt/php/sbt_get_allele_sequence.php")

                        driver.find_element_by_xpath("//select[@id=\'"+gene_index[index]+"\']/option[@value=\'"+str(allele)+"\']").click() # on choisis l allele pour le gene cible

                        time.sleep(wait) # on attend le resultat
                        dna_seq = driver.find_element_by_name('output_text').get_attribute('value') # on recupere la sequence
                        print gene_index[index],' allele ', allele, ' : ', dna_seq

                        driver.close()

                        wait+=1 # on augmente le temp d attente

                    writef.write(dna_seq+'\n')
                    index+=1

        writef.close()


class MyFrame(Frame):
    """
    Interface permettant a l'utilisateur de choisir un lot et de lancer l execution du programme
    pour la recherche de alleles et ST legio SBT
    """

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.frame = Frame(self, relief=RAISED, borderwidth=1) #le frame container qui contient les widget
        self.frame.pack(fill=BOTH, expand=1)

        self.parent = parent #le root

        self.baseDir=r"S:\partage\BioNumericsAnalyse\LegioSBT\SequencesBrutes\\" # localisation des lots de fichiers multifasta sbt

        #self.baseDir=r"C:\Users\ericfournier\Documents\PYTHONtest\SequencesBrutes\\" # localisation des lots de fichiers multifasta sbt

        self.SelectedDir='' # le lot choisis par l utilisateur

        self.LotDirList=[] # liste des lots de fichiers multifasta sbt

        self.MakeDirList() # faire la liste des lots de fichiers multifasta sbt

        #initialiser les widget
        self.initUI()

    def initUI(self):
        """
        Initialisation des widget dans l'interface self.frame
        """

        self.parent.title("Legionella SBT") # titre de l'interface

        self.style = Style() # style d'interface
        self.style.theme_use("clam")

        self.pack(fill=BOTH, expand=1)

        self.messageLabel = Label(self.frame, text=u'Choisir un num\u00E9ro de lot') # message au dessus du combobox
        self.messageLabel.place(x=28,y=30)

        self.okButton = Button(self, text="Ok",fg='blue',command=self.Execute) # bouton Ok pour lancer le programme
        self.okButton.pack(side=RIGHT, padx=5, pady=5 )

        self.quitButton = Button(self, text="Quit",command=self.quit,fg='blue') # bouton Quit pour quitter le programme
        self.quitButton.pack(side=LEFT, padx=5, pady=5)

        self.Cbvar = StringVar() # stringvar lie au combobox
        self.box = ttk.Combobox(self.frame, textvariable=self.Cbvar, state='readonly') # combobox contenant la liste des lots de fichiers multifasta sbt
        self.box.pack(side=LEFT,padx=30)
        self.box['values'] = sorted(self.LotDirList[0:12],reverse=True)

        try:
            self.box.current(0) # lot par defaut affiche dans le combobox
        except:
            box.showinfo('Information',u'Aucun lot disponible à  traiter')
            exit(0)


        self.pB=ttk.Progressbar(self.frame, orient='horizontal', length=200, mode='determinate', value=0,maximum=100) # la progress bar
        self.pB.place(x=20,y=115)

        self.Pblabel = Label(self.frame, text="Ready") # message au dessus de la progress bar
        self.Pblabel.place(x=18,y=90)

    def Execute(self):
        """
        Executer le programme suite au clic Ok
        """

        response=box.askyesno(title='Blast',message="Voulez-vous vraiment continuer ?",icon='warning') # message de confirmaation

        if response==True:

            self.okButton.config(state=DISABLED) #desactiver le bouton ok apres le clic

            self.Pblabel.config(text='s.v.p attendre...',fg='red') # on met a jour le message de la progress bar
            self.update()

            time.sleep(2) # attendre 2 secondes

            self.Pblabel.config(text='In progess ...',fg='red') # on met a jour le message de la progress bar
            self.update()

            self.SelectedDir=self.box.get() # obtenir le lot choisis par l utilisateur

            self.sbt=SBT(self.SelectedDir) # initilisation d un object SBT

            self.sbt.ReadST() # obtenir les profils ST

            self.sbt.ReadSeq() # lire le fichier multifasta sbt du lot choisis

            self.FindAllele() # obtenir les alleles pour chacune des sequences

            self.ToZeroEmptyAllele() # Mettre a zero les alleles manquantes

            self.sbt.WriteSpecAlleleProfil()  # enregister les ST - alleles pour chacun des specimens

            self.quit() # quitter l interface

            shutil.move(self.baseDir+self.SelectedDir+r"/",self.baseDir+"_"+self.SelectedDir+"/") # ajouter un - devant le nom du repertoire du lot

            box.showinfo('Information',u'Terminé')

        else:
            self.quit()

    def ToZeroEmptyAllele(self):
        '''
        Mettre a zero les alleles manquantes
        '''

        for spec in  self.sbt.spec_seq:
            for gene in self.sbt.SBTgene[0:7]:
                if gene not in self.sbt.spec_allele[spec].keys():
                    self.sbt.spec_allele[spec][gene] = '0'


    def FindAllele(self):
        '''
        On recupere les alleles pour chacun des genes SBT des specimens sur le site www.hpa-bioinformatics.org.uk/legionella
        '''
        br = mechanize.Browser()
        br.set_handle_robots(False)   # no robots
        br.set_handle_refresh(False)

        response = br.open("http://www.hpa-bioinformatics.org.uk/legionella/legionella_sbt/php/display_allele_number_from_sequence.php")
        #print response.read()
        br.select_form(name="myForm")

        val=0
        step=float(100/(len(self.sbt.spec_seq)*7))

        for spec in  self.sbt.spec_seq: # pour chacun des specimens du lot
            self.sbt.spec_allele[spec] = {}
            #print "******************* ", spec
            for gene in self.sbt.spec_seq[spec]: # pour chacun des genes sbt
                #print gene
                br.form[gene] = self.sbt.spec_seq[spec][gene] # on met la sequence du gene dans la zone de texte approprie du fomulaire de la page source
                br.submit() # on soumet
                #time.sleep(10)
                #response.read()
                br.select_form(name="myForm")
                #print br.form[gene+'_ST']

                try:
                    self.sbt.spec_allele[spec][re.sub(r'h','',gene)] = br.form[gene+'_ST'] # on recupere l allele resultant. On supprime le h de neuAh.
                except: # probleme avec l une des sequences; aucun allele trouve du a soit; taille incorecte, presence d ambig ou nouvel allele.
                    #box.showinfo('Erreur',u'Problème avec le gène '+gene+u' du spécimen '+ spec)

                    response=box.askyesno(title='Erreur',message=u'Problème avec le gène '+gene+u' du spécimen '+ spec+'.\nVoulez-vous continuer ?',icon='warning') # message de confirmaation

                    if response==True:
                        # on veut eviter d ecraser avec 0 une valeur d allele neuA ou neuAh qui a deja ete enregistree
                        if gene in ['neuAh','neuA']:
                            if 'neuA' in  self.sbt.spec_allele[spec].keys():
                                pass
                            else:
                                self.sbt.spec_allele[spec][re.sub(r'h','',gene)] = '0'

                        else:
                            self.sbt.spec_allele[spec][re.sub(r'h','',gene)] = '0'

                        pass
                    else:
                        self.quit()
                        exit()

                val=val+step # mise a jour de la progress bar
                self.pB.config(value=val)
                self.pB.update()

    def MakeDirList(self):
        """
        Cree la liste des lots de fichiers multifasta a traiter
        """
        self.LotDirList=[]

        for LotDir in os.listdir(self.baseDir):
            if re.search('^\d{8}',LotDir):
                self.LotDirList.append(LotDir)


#le frame racine
root = Tk()
#ajustement des dimensions de la fenetre
root.geometry("450x180+600+600")
#le frame mere
app = MyFrame(root)
#ouvir la fenetre resultante
root.mainloop()


