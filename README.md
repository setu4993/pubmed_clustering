# PubMed Clustering

A project to cluster documents from PubMed based on the disease concepts described in them. The input to this class is a a list of PMIDs, and the output is the clusters for each of the PubMed documents in the input.

The steps followed in this project are:
1. Fetch XML documents from PubMed
2. Parse XML documents to extract title and abstract
3. Pre-process the fetched documents to extract disease words
4. Create a document matrix based on TF-IDF
5. Perform clustering
6. Evaluate results

### Requirements

Environment:
- [Python 3.7](https://www.python.org/downloads/source/) or [Anaconda](https://www.anaconda.com/download/#macos) (preferred) for Python 3.7.

Python libraries (with `pip`):
- [BioPython](http://biopython.org/DIST/docs/install/Installation.html) (for fetching PubMed data): `pip install -U biopython`
- [NumPy](https://scipy.org/install.html) (for matrix operations): `pip install -U numpy scipy matplotlib`
- [SciKit-Learn](https://scikit-learn.org/stable/install.html) (for clustering algorithms): `pip install -U scikit-learn`
- [Gensim](https://radimrehurek.com/gensim/install.html) (for TF-IDF model): `pip install -U gensim`
- [BeautifulSoup4](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup) (to parse XML): `pip install -U beautifulsoup4`

Other tools:
- [UMLS MetaMap](https://metamap.nlm.nih.gov/Installation.shtml) (to extract disease words)
*Note:* A Free UMLS account is required to download the MetaMap binary from the installation link.


### Pipeline
1. The first step is to fetch the documents from PubMed and parse them. They are parsed using BeautifulSoup's XML parser, and the text in the `ArticleTitle` and `Abstract` tags is extracted.
2. The extracted text is then piped to MetaMap. MetaMap queries are run only to identify words from the articles that are disease words. For this, the MetaMap query is limited to show results from `disease or syndrome` and `neoplatic process` semantic types. At the end, the ontology the disease concept is matched to is used instead of the original word that appears in the article. This helps in standardizing concepts and replacing acronyms with the concept.
3. TF-IDF is used to generate the weights for each word in each document. The words appearing in fewer than 3 documents are dropped since they produce noise. Using the TF-IDF weights, the document matrix is generated which is used as an input to the clustering algorithm.
4. The clustering algorithm used is [Affinity Propagation](https://scikit-learn.org/stable/modules/clustering.html#affinity-propagation). The damping factor for Affinity Propagation was set to minimum (`0.5`), after experimenting with various values and observing results and the iterations were set to `500` to ensure it runs to completion and is not terminated. The reasons for choosing Affinity Propagation over other clustering algorithms were the lack of availability of number of clusters, small number of samples and fast convergence. However, this algorithm is not scalable for a large number of samples. If the number of clusters is provided or internal or external metrics are used for evaluation, k-means or hierarchical clustering are better alternatives.
5. The evaluation of the clustering, provided labels, is performed using [purity](https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html) and [F-measure](https://en.wikipedia.org/wiki/F1_score).


### Execution
An example script is provided as a part of `main.py`. This script contains the sample parameters to be passed to the class and function. An example of how to use the provided project is also shown in the Jupyter Notebook `example.ipynb`.

The class `PubMedClustering` contains the methods and variables that are used in executing the whole pipeline. During initialization, the following parameters are to be passed to the function:
- `pubmed_ids` (required): This parameter contains either a list of PMID strings, a string containing comma-separated PMIDs or the disk location to a text file containing a list of PMIDs. If it is a disk location, ensure the `is_file` flag is set to `True`.
- `is_file` (optional): Boolean flag to specify if `pubmed_ids` is a text file that should be read. Default: `False`.
- `metamap` (optional): Location of MetaMap binary. Without MetaMap, the pre-processing will not run and the program will terminate. Default: `/opt/public_mm/bin/metamap16`.
- `email` (optional): Email ID. Required to query PubMed using the BioPython library. Default: `"Your.Name.Here@example.org"`.
- `labels` (optional): Either a dictionary containing PMID -> cluster pairs or a location to a text file that of the form `PMID\tlabel`. If it is a text file, ensure the flag `labels_is_file` is set to `True`. Default: `None`.
- `labels_is_file` (optional): Boolean flag to specify if `pubmed_ids` is a text file that should be read. Default: `False`.


After initializing the class, use the `run()` method of the class to execute the whole pipeline. Upon completion, the results can be seen by querying the `clustering_results_dict` property (a dictionary of PMID -> cluster). If the labels are provided, the evaluation results can be seen by querying the `purity`, `precision`, `recall` and `f_measure` properties of the object.

### Notes
- *Running time:* The running time of the pipeline is `O(N^2 T)`, where `N` is the number of samples and `T` is the iterations because of the clustering algorithm used. All of the other steps are `O(N)`, with the exception of MetaMap queries. Running MetaMap queries depends on the hyper-parameters chosen in MetaMap's execution, but in general, it is slow to execute.
- *Performance Evaluation:* Purity and F-measure (both implemented) are metrics chosen to evaluate the clustering performance.
- *Larger datasets:* Given a larger dataset, a [combination of TF-IDF with word embeddings](https://applied-informatics-j.springeropen.com/articles/10.1186/s40535-018-0055-8) would improve the clustering results. The requirement for a larger dataset is necessary to train a word embeddings model and capture syntactic variations in the dataset.
