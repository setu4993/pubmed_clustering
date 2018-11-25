import os
from pubmed_clustering import fetch_pubmed_documents
import logging
logger = logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s')


if __name__ == '__main__':
    documents = fetch_pubmed_documents(os.path.join('./data', 'pmids_gold_set_unlabeled.txt'), file=True)
    print(documents)
