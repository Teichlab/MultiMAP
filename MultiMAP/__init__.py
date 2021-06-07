import scipy
import numpy as np
from MultiMAP.matrix import MultiMAP, tfidf
#you don't need these if going for MultiMAP.matrix functions
try:
	import anndata
except ImportError:
	pass
try:
	import scanpy as sc
except ImportError:
	pass

def TFIDF_LSI(adata, n_comps=50, binarize=True, random_state=0):
	'''
	Computes LSI based on a TF-IDF transformation of the data. Putative dimensionality 
	reduction for scATAC-seq data prior to MultiMAP. Adds an ``.obsm['X_lsi']`` field to 
	the object it was ran on.
	
	Input
	-----
	adata : ``AnnData``
		The object to run TFIDF + LSI on. Will use ``.X`` as the input data.
	n_comps : ``int``
		The number of components to generate. Default: 50
	binarize : ``bool``
		Whether to binarize the data prior to the computation. Often done during scATAC-seq 
		processing. Default: True
	random_state : ``int``
		The seed to use for randon number generation. Default: 0
	'''
	
	#this is just a very basic wrapper for the non-adata function
	if scipy.sparse.issparse(adata.X):
		adata.obsm['X_lsi'] = tfidf(adata.X.todense(), n_components=n_comps, binarize=binarize, random_state=random_state)
	else:
		adata.obsm['X_lsi'] = tfidf(adata.X, n_components=n_comps, binarize=binarize, random_state=random_state)

def Wrapper(flagged, use_reps, embedding, **kwargs):
	'''
	A function that computes the paired PCAs between the datasets to integrate, calls MultiMAP
	proper, and returns a  (parameters, connectivities, embedding) tuple. Embedding optional
	depending on ``embedding``.
	
	Input
	-----
	flagged : list of ``AnnData``
		Preprocessed objects to integrate. Need to have the single-dataset DRs computed at 
		this stage. Need to have ``.obs[\'multimap_index\']`` defined, incrementing integers
		matching the object's index in the list. Both ``Integrate()`` and ``Batch()`` make 
		these.
	
	All other arguments as described in ``MultiMAP.Integration()``.
	'''
	#MultiMAP wants the shared PCAs delivered as a dictionary, with the subset indices 
	#tupled up as a key. let's make that then
	joint = {}
	#process all dataset pairs
	for ind1 in np.arange(len(flagged)-1):
		for ind2 in np.arange(ind1+1, len(flagged)):
			subset = (ind1, ind2)
			#collapse into a single object and run a PCA
			adata = flagged[ind1].concatenate(flagged[ind2], join='inner')
			sc.tl.pca(adata)
			#preserve space by deleting the intermediate object and just keeping its PCA
			#and multimap index thing
			X_pca = adata.obsm['X_pca'].copy()
			multimap_index = adata.obs['multimap_index'].values
			del adata
			#store the results in joint, which involves some further acrobatics
			joint[subset] = []
			#extract the coordinates for this particular element in the original list, using 
			#the multimap_index .obs column we created before. handy!
			for i in subset:
				joint[subset].append(X_pca[multimap_index == i, :])
	
	#with the joint prepped, we just need to extract the primary dimensionality reductions 
	#and we're good to go here
	Xs = []
	for adata, use_rep in zip(flagged, use_reps):
		Xs.append(adata.obsm[use_rep])
	
	#and with that, we're now truly free to call the MultiMAP function
	#need to negate embedding and provide that as graph_only for the function to understand
	mmp = MultiMAP(Xs=Xs, joint=joint, graph_only=(not embedding), **kwargs)
	
	#and that's it. spit this out for the other wrappers to use however
	return mmp

def Integration(adatas, use_reps, scale=True, embedding=True, **kwargs):
	'''
	Run MultiMAP to integrate a number of AnnData objects from various multi-omics experiments
	into a single joint dimensionally reduced space. Returns a joint object with the resulting 
	embedding stored in ``.obsm[\'X_multimap\']`` (if instructed) and appropriate graphs in 
	``.obsp``. The final object will be a concatenation of the individual ones provided on 
	input, so in the interest of ease of exploration it is recommended to have non-scaled data 
	in ``.X``.
	
	Input
	-----
	adatas : list of ``AnnData``
		The objects to integrate. The ``.var`` spaces will be intersected across subsets of 
		the objects to compute shared PCAs, so make sure that you have ample features in 
		common between the objects. ``.X`` data will be used for computation.
	use_reps : list of ``str``
		The ``.obsm`` fields for each of the corresponding ``adatas`` to use as the 
		dimensionality reduction to represent the full feature space of the object. Needs 
		to be precomputed and present in the object at the time of calling the function.
	scale : ``bool``, optional (default: ``True``)
		Whether to scale the data to N(0,1) on a per-dataset basis prior to computing the 
		cross-dataset PCAs. Improves integration.
	embedding : ``bool``, optional (default: ``True``)
		Whether to compute the MultiMAP embedding. If ``False``, will just return the graph,
		which can be used to compute a regular UMAP. This can produce a manifold quicker,
		but at the cost of accuracy.
	n_neighbors : ``int`` or ``None``, optional (default: ``None``)
		The number of neighbours for each node (data point) in the MultiGraph. If ``None``, 
		defaults to 15 times the number of input datasets.
	n_components : ``int`` (default: 2)
		The number of dimensions of the MultiMAP embedding.
	strengths: ``list`` of ``float`` or ``None`` (default: ``None``)
		The relative contribution of each dataset to the layout of the embedding. The 
		higher the strength the higher the weighting of its cross entropy in the layout loss. 
		If provided, needs to be a list with one 0-1 value per dataset; if ``None``, defaults 
		to 0.5 for each dataset.
	cardinality : ``float`` or ``None``, optional (default: ``None``)
		The target sum of the connectivities of each neighbourhood in the MultiGraph. If 
		``None``, defaults to ``log2(n_neighbors)``.
	
	The following parameter definitions are sourced from UMAP 0.5.1:
	
	n_epochs : int (optional, default None)
		The number of training epochs to be used in optimizing the
		low dimensional embedding. Larger values result in more accurate
		embeddings. If None is specified a value will be selected based on
		the size of the input dataset (200 for large datasets, 500 for small).
	init : string (optional, default 'spectral')
		How to initialize the low dimensional embedding. Options are:
			* 'spectral': use a spectral embedding of the fuzzy 1-skeleton
			* 'random': assign initial embedding positions at random.
			* A numpy array of initial embedding positions.
	min_dist : float (optional, default 0.1)
		The effective minimum distance between embedded points. Smaller values
		will result in a more clustered/clumped embedding where nearby points
		on the manifold are drawn closer together, while larger values will
		result on a more even dispersal of points. The value should be set
		relative to the ``spread`` value, which determines the scale at which
		embedded points will be spread out.
	spread : float (optional, default 1.0)
		The effective scale of embedded points. In combination with ``min_dist``
		this determines how clustered/clumped the embedded points are.
	set_op_mix_ratio : float (optional, default 1.0)
		Interpolate between (fuzzy) union and intersection as the set operation
		used to combine local fuzzy simplicial sets to obtain a global fuzzy
		simplicial sets. Both fuzzy set operations use the product t-norm.
		The value of this parameter should be between 0.0 and 1.0; a value of
		1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy
		intersection.
	local_connectivity : int (optional, default 1)
		The local connectivity required -- i.e. the number of nearest
		neighbors that should be assumed to be connected at a local level.
		The higher this value the more connected the manifold becomes
		locally. In practice this should be not more than the local intrinsic
		dimension of the manifold.
	a : float (optional, default None)
		More specific parameters controlling the embedding. If None these
		values are set automatically as determined by ``min_dist`` and
		``spread``.
	b : float (optional, default None)
		More specific parameters controlling the embedding. If None these
		values are set automatically as determined by ``min_dist`` and
		``spread``.
	'''
	
	#the main thing will be pulling out the various subsets of the adatas, sticking them 
	#together, running joint PCAs, and then splitting up the joint PCAs into datasets of 
	#origin. to do so, let's introduce a helper .obs column in copied versions of adatas
	flagged = []
	for i, adata in enumerate(adatas):
		flagged.append(adata.copy())
		#while we're at it, may as well potentially scale our data copy
		if scale:
			sc.pp.scale(flagged[-1])
		flagged[-1].obs['multimap_index'] = i
	
	#call the wrapper. returns (params, connectivities, embedding), with embedding optional
	mmp = Wrapper(flagged=flagged, use_reps=use_reps, embedding=embedding, **kwargs)
	
	#make one happy collapsed object and shove the stuff in correct places
	#outer join to capture as much gene information as possible for annotation
	adata = anndata.concat(adatas, join='outer')
	if embedding:
		adata.obsm['X_multimap'] = mmp[2]
	#the graph is weighted, the higher the better, 1 best. sounds similar to connectivities
	#TODO: slot distances into .obsp['distances']
	adata.obsp['connectivities'] = mmp[1]
	#set up .uns['neighbors'], setting method to umap as these are connectivities
	adata.uns['neighbors'] = {}
	adata.uns['neighbors']['params'] = mmp[0]
	adata.uns['neighbors']['params']['method'] = 'umap'
	adata.uns['neighbors']['distances_key'] = 'distances'
	adata.uns['neighbors']['connectivities_key'] = 'connectivities'
	return adata

def Batch(adata, batch_key='batch', scale=True, embedding=True, dimred_func=None, rep_name='X_pca', **kwargs):
	'''
	Run MultiMAP to correct batch effect within a single AnnData object. Loses the flexibility 
	of individualised dimensionality reduction choices, but doesn't require a list of separate 
	objects for each batch/dataset to integrate. Runs PCA on a per-batch/dataset basis prior 
	to performing an analysis analogous to  ``Integration()``. Adds appropriate ``.obsp`` graphs 
	and ``.obsm[\'X_multimap\']`` (if instructed) to the input.
	
	Input
	-----
	adata : ``AnnData``
		The object to process. ``.X`` data will be used in the computation.
	batch_key : ``str``, optional (default: "batch")
		The ``.obs`` column of the input object with the categorical variable defining the 
		batch/dataset grouping to integrate on.
	scale : ``bool``, optional (default: ``True``)
		Whether to scale the data to N(0,1) on a per-dataset basis prior to computing the 
		cross-dataset PCAs. Improves integration.
	embedding : ``bool``, optional (default: ``True``)
		Whether to compute the MultiMAP embedding. If ``False``, will just return the graph,
		which can be used to compute a regular UMAP. This can produce a manifold quicker,
		but at the cost of accuracy.
	dimred_func : function or ``None``, optional (default: ``None``)
		The function to use to compute dimensionality reduction on a per-dataset basis. Must 
		accept an ``AnnData`` on input and modify it by inserting its dimensionality reduction 
		into ``.obsm``. If ``None``, ``scanpy.tl.pca()`` will be used.
	rep_name : ``str``, optional (default: "X_pca")
		The ``.obsm`` field that the dimensionality reduction function stores its output under.
	
	All other arguments as described in ``Integration()``.
	'''
	
	#as promised in the docstring, set dimred_func to scanpy PCA if not provided
	if dimred_func is None:
		dimred_func = sc.tl.pca
	
	#essentially what this function does is preps data to run through the other wrapper
	#so what needs to happen is the object needs to be partitioned up, have DR ran,
	#and passed as a list to the wrapper function
	flagged = []
	use_reps = []
	for i,batch in enumerate(np.unique(adata.obs[batch_key])):
		#extract the single batch data
		flagged.append(adata[adata.obs[batch_key]==batch].copy())
		#potentially scale
		if scale:
			sc.pp.scale(flagged[-1])
		#and run DR
		dimred_func(flagged[-1])
		#and stick on the index for multimap to pull stuff apart later
		flagged[-1].obs['multimap_index'] = i
		#and add an entry to the list of .obsm keys for the other function
		use_reps.append(rep_name)
	
	#call the wrapper. returns (params, connectivities, embedding), with embedding optional
	mmp = Wrapper(flagged=flagged, use_reps=use_reps, embedding=embedding, **kwargs)
	
	#stick stuff where it's supposed to go
	if embedding:
		adata.obsm['X_multimap'] = mmp[2]
	#the graph is weighted, the higher the better, 1 best. sounds similar to connectivities
	#TODO: slot distances into .obsp['distances']
	adata.obsp['connectivities'] = mmp[1]
	#set up .uns['neighbors'], setting method to umap as these are connectivities
	adata.uns['neighbors'] = {}
	adata.uns['neighbors']['params'] = mmp[0]
	adata.uns['neighbors']['params']['method'] = 'umap'
	adata.uns['neighbors']['distances_key'] = 'distances'
	adata.uns['neighbors']['connectivities_key'] = 'connectivities'