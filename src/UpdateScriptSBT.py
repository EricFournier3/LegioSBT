#encoding: iso-8859-15
#//$autorun;event=PostCreateMainWin
import bns,shelve
import os
import os.path
import re
import datetime
import shutil
import time
from xml.dom.minidom import parse

'''
Ce script s'execute automatiquement apres l'ouverture de la bd LegioSBT.
Si necessaire, il met a jour les scripts SBTAssemblyTemplate.py, ImportSpecST.py et UpdateScript.py dans les profils utilisateurs.
Il met egalement a jour les bd de reference(SBT_Ref) et locale(SBT_Locale)
'''

MessageBox=bns.Util.Program.MessageBox 

DatabaseName=bns.Database.Db.Info.Name  # nom de la bd courante
UserName=bns.Database.Db.Info.WindowsLogin  # utilisateur

def UpdateBdUser():
    """
    Mise a jour de SBT_ref et SBT_Locale
    """

    # Formattage de la date de la mise a jour
    dateUpdate=str(datetime.datetime.now())
    dateformat=re.search(r'(\d{4}-\d{2}-\d{2})',dateUpdate)
    dateUpdatePart=dateformat.group(1)
    dateUpdatePart=re.sub(r'-','',dateUpdatePart)

    Entries=bns.Database.Db.Entries # liste des entrees BioNumerics

    bns.Windows.BnsWindow.SetWaitCursor(True) # sablier

    for entry in Entries: # on selectionne exlusivement toutes les sequences de reference
        if re.search(r'^ST_\d+',entry.Key):
            entry.Selection=True
        else:
            entry.Selection=False

    bns.Blast.BlastTools.CreateBlastDb('nucleic acids','SBT_Ref',dateUpdatePart) # mise a jour de SBT_Ref

    for entry in Entries: # on selectionne exclusivement toutes les specimens cliniques
        if entry.Selection == True:
            entry.Selection = False
        else:
            entry.Selection = True

    bns.Blast.BlastTools.CreateBlastDb('nucleic acids','SBT_Locale',dateUpdatePart) # mise a jour de SBT_Locale

    for entry in Entries: # on deselectionne toutes les entrees
        entry.Selection=False


def UpdateScript(ScriptFrom,ScriptTo):
		'''
		Pour la mise a jour des scripts SBTAssemblyTemplate.py, ImportSpecST.py et UpdateScript.py  de S:\partage\BioNumericsAnalyse\Scripts\Python\BioNumerics\LegioSBT\
		vers U:\Citrix\Documents\BioNumerics\Data\LegioSBT\Scripts\Python
		'''

		newScriptFrom=ScriptFrom
		newScriptTo=ScriptTo

		#supprimer les scripts  existants dans newScriptTo
		for script in os.listdir(newScriptTo):
			os.remove(newScriptTo+str(script))

		#copie des nouveaux scripts bacterio de newScriptFrom vers newScriptTo
		for script in os.listdir(newScriptFrom):
			shutil.copy2(newScriptFrom+str(script),newScriptTo)
				
		MessageBox(" " ,'Bienvenue '+UserName+ ' dans la base de données '+DatabaseName+'\n'
		+'Mise à jour des scripts à effectuer. Veuillez svp rouvrir la bd.','information')
		
		# Fermer la fenetre principale; le Script Window en mode debug ou le BioNumerics Main Window qui apparait apres l'ouverture de la bd
		# La fenetre principale a toujours le id 1
		bns.Windows.BnsWindow(1).Close()
				

def CheckForScriptUpdate():
    """
    Verifier si une mise a jour est necessaire
    """

    ScriptFrom = r"S:\partage\BioNumericsAnalyse\Scripts\Python\BioNumerics\legioSBT\\" #  path nouveaux scripts dans BioNumericsAnalyse
    ScriptTo=r"U:\Citrix\Documents\BioNumerics\Data\LegioSBT\Scripts\Python\\"  # path vers les scripts du profil utilisateur

    ListScriptFrom=[]  #  Liste des dates modifs des scripts  dans ScriptFrom
    ListScriptTo=[]  # Liste des dates modifs des scripts  dans ScriptTo

    #faire la liste des dates modif des scripts dans ScriptTo
    for script in os.listdir(ScriptTo):
        ListScriptTo.append(time.ctime(os.path.getmtime(ScriptTo+script)))


    for script in os.listdir(ScriptFrom):
        ListScriptFrom.append(time.ctime(os.path.getmtime(ScriptFrom+script)))

    # si aucune mise a jour des scripts necessaires on peut mettre a jour les bd
    if (ListScriptFrom==ListScriptTo):
        UpdateBdUser()
    else:
        UpdateScript(ScriptFrom,ScriptTo)
		
				
# verifier s'il y a lieu de mettre a jour les scripts	
CheckForScriptUpdate()


