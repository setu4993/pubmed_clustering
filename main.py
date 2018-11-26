import os
from pubmed_clustering import PubMedClustering
import logging

logger = logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s')


if __name__ == '__main__':
    clus = PubMedClustering(os.path.join('./data', 'pmids_gold_set_unlabeled.txt'),
                            metamap='/opt/public_mm/bin/metamap16', email="This.Is.My@Email.org", is_file=True,
                            labels=os.path.join('./data', 'pmids_gold_set_labeled.txt'), labels_is_file=True)
    clus.run()
    print('Clustering results: ', clus.clustering_results)
    print('Purity:', clus.purity)
    print('F-measure:', clus.f_measure)
