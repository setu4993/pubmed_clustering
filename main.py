import os
from pubmed_clustering import PubMedClustering
import logging

logger = logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s')


if __name__ == '__main__':
    clus = PubMedClustering(os.path.join('./data', 'pmids_gold_set_unlabeled.txt'), is_file=True)
    clus.run()
    print(clus.clustering_results)
