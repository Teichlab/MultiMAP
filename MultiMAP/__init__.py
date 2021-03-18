from itertools import chain, combinations
import scipy
import numpy as np
from MultiMAP.matrix import MultiGraph, MultiMAP, tfidf
#you don't need these if going for MultiMAP.matrix functions
try:
	import anndata
except ImportError:
	pass
try:
	import scanpy as sc
except ImportError:
	pass

def powerset(iterable, minlen=0):
	'''
	A function to get all of the subsets of a given set
	
	Input:
		- iterable - the set to find subsets for
		- minlen - the minimum subset size to return
	
	Returns an itertools.chain with tuples of the subsets
	'''
	s = list(iterable)
	return chain.from_iterable(combinations(s, r) for r in range(minlen, len(s)+1))

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

def MultiMAP_Integration(adatas, use_reps, scale=True, **kwargs):
	'''
	Run MultiMAP to integrate a number of AnnData objects from various multi-omics experiments
	into a single joint dimensionally reduced space. Returns a joint object with the resulting 
	embedding stored in ``.obsm[\'X_multimap\']`` and appropriate graphs in ``.obsp``. The 
	final object will be a concatenation of the individual ones provided on input, so in the 
	interest of ease of exploration it is recommended to have non-scaled data in ``.X``.
	
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
	same_genes = True
	for i, adata in enumerate(adatas):
		flagged.append(adata.copy())
		#while we're at it, may as well potentially scale our data copy
		if scale:
			sc.pp.scale(flagged[-1])
		flagged[-1].obs['multimap_index'] = i
	
	#MultiMAP wants the shared PCAs delivered as a dictionary, with the subset indices 
	#tupled up as a key. let's make that then
	joint = {}
	#we just want subsets 2 or bigger
	for subset in powerset(np.arange(len(flagged)), 2):
		#so you can't subset lists by many indices, and a list of AnnDatas doesn't like 
		#the idea of becoming a np.array... engage clumsy resolution!
		subflag = []
		for i in subset:
			subflag.append(flagged[i])
		#collapse into a single object and run a PCA
		adata = anndata.concat(subflag, join='inner')
		sc.pp.pca(adata)
		#store the results in joint, which involves some further acrobatics
		joint[subset] = []
		#extract the coordinates for this particular element in the original list, using 
		#the multimap_index .obs column we created before. handy!
		for i in subset:
			asub = adata[adata.obs['multimap_index'] == i]
			joint[subset].append(asub.obsm['X_pca'])
	
	#with the joint prepped, we just need to extract the primary dimensionality reductions 
	#and we're good to go here
	Xs = []
	for adata, use_rep in zip(adatas, use_reps):
		Xs.append(adata.obsm[use_rep])
	
	#and with that, we're now truly free to call the MultiMAP functions
	#TODO: catch another argument once Mika distances appear
	embed, connectivities, params = MultiMAP(Xs=Xs, joint=joint, **kwargs)
	
	#make one happy collapsed object and shove the stuff in correct places
	#outer join to capture as much gene information as possible for annotation
	adata = anndata.concat(adatas, join='outer')
	adata.obsm['X_multimap'] = embed
	#the graph is weighted, the higher the better, 1 best. sounds similar to connectivities
	#TODO: slot distances into .obsp['distances']
	adata.obsp['connectivities'] = connectivities
	#set up .uns['neighbors'], setting method to umap as these are connectivities
	adata.uns['neighbors'] = {}
	adata.uns['neighbors']['params'] = params
	adata.uns['neighbors']['params']['method'] = 'umap'
	adata.uns['neighbors']['distances_key'] = 'distances'
	adata.uns['neighbors']['connectivities_key'] = 'connectivities'
	return adata

def MultiMAP_Batch(adata, batch_key='batch', scale=True, dimred_func=None, rep_name='X_pca', **kwargs):
	'''
	Run MultiMAP to correct batch effect within a single AnnData object. Loses the flexibility 
	of individualised dimensionality reduction choices, but doesn't require a list of separate 
	objects for each batch/dataset to integrate. Runs PCA on a per-batch/dataset basis prior 
	to calling ``MultiMAP_Integration()`` internally. Adds ``.obsm['X_multimap']`` and 
	appropriate ``.obsp`` graphs to the input.
	
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
	dimred_func : function or ``None``, optional (default: ``None``)
		The function to use to compute dimensionality reduction on a per-dataset basis. Must 
		accept an ``AnnData`` on input and modify it by inserting its dimensionality reduction 
		into ``.obsm``. If ``None``, ``scanpy.tl.pca()`` will be used.
	rep_name : ``str``, optional (default: "X_pca")
		The ``.obsm`` field that the dimensionality reduction function stores its output under.
	
	All other arguments as described in ``MultiMAP_Integration()``.
	'''
	
	#as promised in the docstring, set dimred_func to scanpy PCA if not provided
	if dimred_func is None:
		dimred_func = sc.tl.pca
	
	#essentially what this function does is preps data to run through the other wrapper
	#so what needs to happen is the object needs to be partitioned up, have DR ran,
	#and passed as a list to the integration function
	adatas = []
	use_reps = []
	for batch in np.unique(adata.obs[batch_key]):
		#extract the single batch data
		adatas.append(adata[adata.obs[batch_key]==batch].copy())
		#potentially scale
		if scale:
			sc.pp.scale(adatas[-1])
		#and run DR
		dimred_func(adatas[-1])
		#and add an entry to the list of .obsm keys for the other function
		use_reps.append(rep_name)
	
	#and that's it, really. call the other function
	#set scale to False - if this function's scale was True, the data's already scaled
	#if it was set to False, then we're not scaling anyway
	mmp = MultiMAP_Integration(adatas=adatas, use_reps=use_reps, scale=False, **kwargs)
	
	#this made us a new object, but we just want to copy stuff into the original one
	adata.obsm['X_multimap'] = mmp.obsm['X_multimap']
	#TODO: catch distances once those are generated
	adata.obsp['connectivities'] = mmp.obsp['connectivities']
	adata.uns['neighbors'] = mmp.uns['neighbors']