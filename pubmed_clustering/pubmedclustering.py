from .pubmed_docs_fetch import fetch_pubmed_documents
from .metamap_operations import documents_disease_list
import numpy as np
from gensim.models import TfidfModel
from gensim.corpora import Dictionary
import logging
from sklearn.cluster import AffinityPropagation

logger = logging.getLogger()


class PubMedClustering:
    def __init__(self, pubmed_ids, metamap='/opt/public_mm/bin/metamap16', email="Your.Name.Here@example.org",
                 is_file=False, labels=None, labels_is_file=None):
        self.pubmed_ids = pubmed_ids
        self.metamap_location = metamap
        self.email = email
        self.is_file = is_file
        self.target_labels = labels
        self.labels_is_file = labels_is_file

        self.documents = None
        self.pmids = None

        self.disease_words = None
        self.metamap_json = None

        self.dictionary = None
        self.corpus = None
        self.tfidf_model = None

        self.affinity_propagation = None
        self.clustering_results = None
        self.clustering_results_dict = None
        self.cluster_wise_results = None

        self.ground_truth = None
        self.categories = None

        self.purity = None
        self.precision = None
        self.recall = None
        self.f_measure = None

    def run(self):
        self.documents, self.pmids = fetch_pubmed_documents(self.pubmed_ids, email=self.email, is_file=self.is_file)
        self.disease_words, self.metamap_json = documents_disease_list(self.documents)
        self.__create_tfidf_model()
        self.__calc_document_matrix()
        self.__clustering()
        if self.target_labels:
            self.__original_labels(self.target_labels, is_file=self.labels_is_file)
            self.__purity()
            self.__f_measure()

    def __create_tfidf_model(self):
        logging.info('Creating dictionary')
        self.dictionary = Dictionary(self.disease_words)
        logging.info('Filtering low occurence words from dictionary')
        self.dictionary.filter_extremes(no_below=3, no_above=1.0)
        self.dictionary.compactify()
        logging.info('Creating a corpus')
        self.corpus = [self.dictionary.doc2bow(line) for line in self.disease_words]
        logging.info('Creating a TF-IDF model')
        self.tfidf_model = TfidfModel(self.corpus)

    def __calc_document_matrix(self):
        logging.info('Calculating document matrix')

        self.document_matrix = np.zeros((len(self.disease_words), len(self.dictionary)))
        for i, disease_words in enumerate(self.corpus):
            if i % 10 == 0:
                logging.info('Document vector for %d of %d documents calculated', i, len(self.documents))
            tfidf_score = self.tfidf_model[disease_words]
            for (dct_id, score) in tfidf_score:
                self.document_matrix[i][dct_id] = score
        logging.info('Document matrix calculated')

    def __clustering(self):
        logging.info('Performing clustering on the document matrix')
        self.affinity_propagation = AffinityPropagation(damping=0.5, max_iter=500).fit(self.document_matrix)
        logging.info('Clustering performed')
        self.clustering_results = self.affinity_propagation.labels_
        self.clustering_results_dict = {pmid: label for pmid, label in zip(self.pmids, self.clustering_results)}
        self.total_clusters = len(self.affinity_propagation.cluster_centers_indices_)
        logging.info('Clustering results saved to object')

    def __original_labels(self, ground_truth, is_file=False):
        logging.info('Reading ground truth labels')
        if is_file:
            doc_lines = open(ground_truth, 'r').readlines()
            ground_truth = self.__read_labels(doc_lines)

        if type(ground_truth) is dict:
            self.ground_truth = ground_truth
            self.categories = set(ground_truth.values())
        else:
            raise TypeError('Incorrect type of `ground_truth` variable')

    def __read_labels(self, doc_lines):
        categories_dict = {}
        pmid_ground_truth = {}
        i = 0
        for line in doc_lines:
            contents = line.strip().split('\t')
            if contents[1] not in categories_dict:
                categories_dict[contents[1]] = i
                i += 1
            pmid_ground_truth[contents[0]] = categories_dict[contents[1]]
        return pmid_ground_truth

    def __purity(self):
        logging.info('Calculating purity')
        self.cluster_wise_results = [{} for _ in range(self.total_clusters)]
        for pmid, cluster in self.clustering_results_dict.items():
            if self.ground_truth[pmid] not in self.cluster_wise_results[cluster]:
                self.cluster_wise_results[cluster][self.ground_truth[pmid]] = 0
            self.cluster_wise_results[cluster][self.ground_truth[pmid]] += 1

        cluster_max = [max(cluster_result.values()) for cluster_result in self.cluster_wise_results]
        self.purity = sum(cluster_max) / len(self.clustering_results_dict)
        logging.info('Purity calculated and stored')

    def __f_measure(self):
        logging.info('Calculating F-measure')

        self.precision = np.mean([max(cluster_result.values()) / sum(cluster_result.values()) for cluster_result in
                                  self.cluster_wise_results])
        recall_list = {category: 0 for category in self.categories}
        recall_total = {category: 0 for category in self.categories}
        for cluster_result in self.cluster_wise_results:
            for category, count in cluster_result.items():
                recall_total[category] += count
            max_key = max(cluster_result, key=cluster_result.get)
            recall_list[max_key] += cluster_result[max_key]
        self.recall = np.mean([recall_list[category] / recall_total[category] for category in
                               self.categories])
        self.f_measure = (2 * self.precision * self.recall) / (self.precision + self.recall)
        logging.info('Precision, recall and F-measure calculated and stored')
