##########################################################
#
# Fonctions interactives permettant de faire des choix
#  a l'ecran, en utilisant la souris et des fenetres
#  graphiques qui le font bien !
#
# Creation: 2005/2/16 - A. Sarrat
# Modifie: 2005/2/19 -  F. Schuller
#          .... 2005/10/28
# Last modif: 2005/11/04 - adapted to BoGLi
#

from ppgplot import *
import string


class Fenetre:
  """
  classe Fenetre - parametres et methodes pour les boites et boutons
  attributs:
  int forme    : 0=cercle 1=rectangle 2=rectangle transparent
  list pos     : positions (X,Y) des centres dans la fenetre
  list float/tuple size: rayon ou (largeur,hauteur)
  list label   : messages a apparaitre (vecteur de string) dans ou pres d'une fenetre
  tuple txtpos : position des labels relativement a pos
  int font     : taille des caracteres
  int coltxt   : couleur du texte
  int colfond  : couleur de fond des boutons
  int family   : police de caracteres
  """
  
  def __init__(self,forme=0,pos=[],size=0,label=[],txtpos=(0,0),font=2.,\
  		coltxt=1,colfond=1,family=1):
    self.forme = forme
    self.pos = pos
    self.size = size
    self.label = label
    self.txtpos = txtpos
    self.font = font
    self.coltxt = coltxt
    self.colfond = colfond
    self.family = family

  def dessine(self,new=0):
    """
    method dessine
    
    INP: new : efface la fenetre si non nul
    OUT: aucune
    """
    if new:
      pgeras()
    pgsfs(2)
    # pgrect(0,1000,500,0)
    pgsave()                # sauvegarde les attributs de la figure
    pgbbuf()                # demarre la preparation de la page (batch)
    pgslw(2)                # definit l'epaisseur des lignes
    pgsfs(1)                # definit le style de remplissage
    pgsch(self.font)        # definit la taille des caracteres
    pgscf(self.family)      # definit la fonte

    # Cree une mini-page temporaire de selection
    for i in range(len(self.label)):
      pgsci(self.colfond)       # definit l'index de couleur 1:blanc 2:rouge...
      x0 = self.pos[i][0]
      y0 = self.pos[i][1]
      if self.forme == 1:
        oneSize = self.size[i]
        pgrect(x0-oneSize[0]/2.,x0+oneSize[0]/2., \
            y0-oneSize[1]/2.,y0+oneSize[1]/2.,)
      elif self.forme == 0:
        oneSize = (self.size[i],self.size[i])
        pgcirc(x0,y0,self.size[i])
      pgsci(self.coltxt)
      if '|' in self.label[i]:
        # allow for two words, separated on two lines
        words = string.split(self.label[i],'|')
        pgptxt(x0+self.txtpos[0], y0+self.txtpos[1]-oneSize[1]/4.,\
               0.0, 0.5, words[0])
        pgptxt(x0+self.txtpos[0], y0+self.txtpos[1]+oneSize[1]/4.,\
               0.0, 0.5, words[1])
      else:
        pgptxt(x0+self.txtpos[0], y0+self.txtpos[1], 0.0, 0.5, self.label[i])
    pgebuf()       # fin de la preparation de la page
    pgunsa()       # retourne aux attributs de pgplot


  def saisie(self):
    """
    method saisie
    
    INP: aucune
    OUT: choix : selection (numero du bouton)
    """
    
    # Selection d'un bouton de reponse par l'utilisateur (a la souris)
    again = 1
    xx, yy = 500, 150
    while(again):
      xx,yy,zz = pgcurs(xx,yy)
      choix = -1
      for i in range(len(self.label)):
        x0 = self.pos[i][0]
        y0 = self.pos[i][1]
        if self.forme:
          oneSize = self.size[i]
          if ( abs(xx-x0) <= oneSize[0]/2. and \
               abs(yy-y0) <= oneSize[1]/2.):
            choix = i
            again = 0
        else:
          oneSize = (self.size[i],self.size[i])
          if ((xx-x0)**2 + (yy-y0)**2) <= oneSize**2:
            choix = i
            again = 0
            
    return choix


def fenetreInteractive(x0=600.,x1=675.,y0=400.,y1=480,prevText='',bgCol=0):
  """
  method fenetreInteractive():
  INP: (float) x0,x1,y0,y1: box coordinates
       (str) prevText: previous field value
       (int) bgCol: background color
       
  OUT: (str) value: user input
  """

  box = [x0,x1,y0,y1]
  X = x0 + 10.
  Y = y1 - 10.

  pgsave()                # sauvegarde les attributs de la figure
  pgbbuf()                # demarre la preparation de la page (batch)
  pgsci(16)               # definit l'index de couleur 1:blanc 2:rouge...
  pgslw(2)                # definit l'epaisseur des lignes
  pgsfs(1)                # definit le style de remplissage
  pgsch(2.0)              # definit la taille des caracteres
  pgscf(1)                # definit la fonte

  # black frame
  pgrect(x0,x1,y0,y1)

  # input text in red
  pgsci(2)
  value = pgrstr(X, Y, 0.0, 0.0, prevText, len(prevText), bgCol)
  value = value.rstrip()
  
  pgebuf()       # fin de la preparation de la page
  pgunsa()       # retourne aux attributs de pgplot

  return value


def pgrstr(X, Y, ANGLE, FJUST, TEXT, LSTR, BCI):
  """
  method pgrstr(X, Y, ANGLE, FJUST, TEXT, LSTR, BCI)
  Lit une chaine de caractere dans la fenetre courante
  INP: X, Y, ANGLE, FJUST = position angle et justification du texte
       BCI = couleur de la fenetre
  INP/OUT: TEXT, LSTR = texte et longueur de la chaine TEXT
  """

  CI = pgqci()
  again = 1
  while(again):
    # Draw current string
    if ( LSTR > 0 ):
        pgptxt(X, Y, ANGLE, FJUST, TEXT[0:LSTR])
        JUNK = pgqtxt(X, Y, ANGLE, FJUST, TEXT[0:LSTR])
        XCUR = JUNK[2][0]
        YCUR = JUNK[2][1]
    else:
        XCUR = X
        YCUR = Y

    # Read a character
    XX,YY,CH = pgband(0, 1, XCUR, YCUR)

    # Erase old string
    pgsci(BCI)
    if (LSTR > 0):
        pgptxt(X, Y, ANGLE, FJUST, TEXT[0:LSTR])
    pgsci(CI)

    # Avoid problem with PGPLOT escape character
    if (CH == chr(92)):
        CH = '*'
        
    # Backspace (ctrl H) or delete remove entire string
    if ( ord(CH) == 8 or ord(CH) == 127 ):
        if ( LSTR > 0 ):
            LSTR = 0
            TEXT = ""

    # Any other non-printing character terminates input
    elif ( ord(CH) < 32 ):
        if ( LSTR > 0 ):
            pgptxt(X, Y, ANGLE, FJUST, TEXT[0:LSTR])
        again = 0

    # Otherwise, add character to string if there is room
    else:
        LSTR = LSTR + 1
        TEXT = TEXT + CH

  return TEXT

## @}
