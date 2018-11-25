from Bio import Entrez
from bs4 import BeautifulSoup
import logging

logger = logging.getLogger()


def fetch_pubmed_documents(pmids, email="Your.Name.Here@example.org", file=False):
    if file:
        pmids = open(pmids, 'r').readlines()

    if type(pmids) is list:
        pmids = ",".join(pmids)

    return _parse_pubmed_xml(_query_pubmed(pmids, email))


def _query_pubmed(pmids, email="Your.Name.Here@example.org"):
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract")
    return handle.read()


def _parse_pubmed_xml(pubmed_xml):
    soup = BeautifulSoup(pubmed_xml, 'xml')

    documents = []

    for i, record in enumerate(soup.find_all('PubmedArticle')):
        try:
            pmid = record.MedlineCitation.PMID.text
        except AttributeError:
            pmid = None
        try:
            title = record.MedlineCitation.Article.ArticleTitle.text
        except AttributeError:
            title = ''
        try:
            abstract = record.MedlineCitation.Article.Abstract.text
        except AttributeError:
            abstract = ''
        documents.append({'pmid': pmid, 'title': title, 'abstract': abstract, 'all_text': title + ' ' + abstract})

    return documents
