from .pubmed_docs_fetch import fetch_pubmed_documents
from .metamap_operations import documents_disease_list
import numpy as np
from gensim.models import TfidfModel
from gensim.corpora import Dictionary
import logging
from sklearn.cluster import AffinityPropagation

logger = logging.getLogger()


class PubMedClustering:
    def __init__(self, pubmed_ids, metamap='/opt/public_mm/bin/metamap16', email="Your.Name.Here@example.org", is_file=False):
        self.pubmed_ids = pubmed_ids
        self.metamap_location = metamap
        self.email = email
        self.is_file = is_file

        self.documents = None
        self.pmids = None

        self.disease_words = None
        self.metamap_json = None

        self.dictionary = None
        self.corpus = None
        self.tfidf_model = None

        self.affinity_propagation = None
        self.clustering_results = None

    def run(self):
        self.documents, self.pmids = fetch_pubmed_documents(self.pubmed_ids, email=self.email, is_file=self.is_file)
        self.disease_words, self.metamap_json = documents_disease_list(self.documents)
        self._create_tfidf_model()
        self._calc_document_matrix()
        self._clustering()

    def _create_tfidf_model(self):
        logging.info('Creating dictionary')
        self.dictionary = Dictionary(self.disease_words)
        logging.info('Filtering low occurence words from dictionary')
        self.dictionary.filter_extremes(no_below=3, no_above=1.0)
        self.dictionary.compactify()
        logging.info('Creating a corpus')
        self.corpus = [self.dictionary.doc2bow(line) for line in self.disease_words]
        logging.info('Creating a TF-IDF model')
        self.tfidf_model = TfidfModel(self.corpus)

    def _calc_document_matrix(self):
        logging.info('Calculating document matrix')

        self.document_matrix = np.zeros((len(self.disease_words), len(self.dictionary)))
        for i, disease_words in enumerate(self.corpus):
            if i % 10 == 0:
                logging.info('%d of %d documents processed with MetaMap', i, len(self.documents))
            tfidf_score = self.tfidf_model[disease_words]
            for (dct_id, score) in tfidf_score:
                self.document_matrix[i][dct_id] = score
        logging.info('Document matrix calculated')

    def _clustering(self):
        logging.info('Performing clustering on the document matrix')
        self.affinity_propagation = AffinityPropagation(damping=0.5, max_iter=500).fit(self.document_matrix)
        logging.info('Clustering performed')
        self.clustering_results = self.affinity_propagation.labels_
        logging.info('Clustering results saved to class')
