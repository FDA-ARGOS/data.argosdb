'''
This is a model script needed to search and retrieve a list of PMIDs relating to a user-specified gene. 
'''

from Bio import Entrez

Entrez.email = 'email@email.com'
handle = Entrez.esearch(db='pubmed', term=('Salmonella enterica[title]'+ 'substitution[title]' + 'ramR'))
record = Entrez.read(handle)
print("Salmonella enterica PMIDs:\n",record['IdList'])
