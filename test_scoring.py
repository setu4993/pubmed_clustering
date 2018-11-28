from pubmed_clustering import PubMedClustering
import unittest


class TestScoringMethods(unittest.TestCase):

    def create_pubmed_test_object(pmids, labels, clusters):
        clus = PubMedClustering(pmids, labels=labels)
        clus._PubMedClustering__original_labels(clus.target_labels, clus.labels_is_file)
        clus.clustering_results_dict = clusters
        clus.total_clusters = len(set(clusters.values()))
        clus._PubMedClustering__purity()
        clus._PubMedClustering__f_measure()
        return clus

    def test_single_cluster_all_correct(self):
        pmids = ['3242', '324234', '4356', '2355', '32626', '35646']
        pmids_labels = {'3242': 0, '324234': 0, '4356': 0, '2355': 0, '32626': 0, '35646': 0}
        pmids_clusters = pmids_labels
        clus = self.create_pubmed_test_object(pmids, pmids_labels, pmids_clusters)
        self.assertAlmostEqual(clus.purity, 1.0)
        self.assertAlmostEqual(clus.f_measure, 1.0)
        self.assertAlmostEqual(clus.recall, 1.0)
        self.assertAlmostEqual(clus.precision, 1.0)

    def test_two_clusters_all_correct(self):
        pmids = ['3242', '324234', '4356', '2355', '32626', '35646']
        pmids_labels = {'3242': 0, '324234': 0, '4356': 0, '2355': 0, '32626': 0, '35646': 0}
        pmids_clusters = {'3242': 0, '324234': 0, '4356': 0, '2355': 1, '32626': 1, '35646': 1}
        clus = self.create_pubmed_test_object(pmids, pmids_labels, pmids_clusters)
        self.assertAlmostEqual(clus.purity, 1.0)
        self.assertAlmostEqual(clus.f_measure, 1.0)
        self.assertAlmostEqual(clus.recall, 1.0)
        self.assertAlmostEqual(clus.precision, 1.0)

    def test_two_clusters_one_incorrect(self):
        pmids = ['3242', '324234', '2355', '32626', '35646']
        pmids_labels = {'3242': 0, '324234': 0, '2355': 1, '32626': 0, '35646': 0}
        pmids_clusters = {'3242': 1, '324234': 1, '2355': 0, '32626': 0, '35646': 0}
        clus = self.create_pubmed_test_object(pmids, pmids_labels, pmids_clusters)
        self.assertAlmostEqual(clus.purity, 0.8)
        self.assertAlmostEqual(clus.f_measure, 0.625)
        self.assertAlmostEqual(clus.recall, 0.5)
        self.assertAlmostEqual(clus.precision, 0.833333333333)

    def test_two_clusters_two_two_one_incorrect(self):
        pmids = ['3242', '324234', '2355', '32626', '35646']
        pmids_labels = {'3242': 0, '324234': 0, '2355': 1, '32626': 1, '35646': 1}
        pmids_clusters = {'3242': 1, '324234': 1, '2355': 1, '32626': 0, '35646': 0}
        clus = self.create_pubmed_test_object(pmids, pmids_labels, pmids_clusters)
        self.assertAlmostEqual(clus.purity, 0.8)
        self.assertAlmostEqual(clus.f_measure, 0.833333333333)
        self.assertAlmostEqual(clus.recall, 0.833333333333)
        self.assertAlmostEqual(clus.precision, 0.833333333333)


if __name__ == '__main__':
    unittest.main()
