from sklearn.decomposition import PCA
from sklearn.manifold import TSNE;


data; #Your Data Matrix - rows are Genes, columns are Samples

#PCA
#Documentation here - http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
model = PCA(n_components=2);
result_pca = model.fit_transform(data.T);  #Note the transpose
#Now result_pca is of shape NUM_SAMPLES x 2

#tSNE
#Documentation here - http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
model = TSNE(n_components=2, perplexity=30.0, metric="euclidean", learning_rate = 100, early_exaggeration=4.0);
result_tsne = model.fit_transform(data.T);  #Note the transpose

#Now result_tsne is of shape NUM_SAMPLES x 2
