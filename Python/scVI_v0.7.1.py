import argparse
import os
import pandas as pd
import numpy as np
#import torch
import warnings; warnings.simplefilter('ignore')
import scanpy as sc
import re
import sys
import scvi


'''
os.chdir changes the working directory so that the code should work the same as when I run it from command line.
This works as I keep my scripts in the ./scripts directory from the projects home directory ./ 
'''
# os.chdir(os.getcwd() + '/scripts') 
# os.chdir(os.getcwd() + '/github/Python') 
# os.chdir('/local/data1/user/sanli71/CIA_project/github/CIA_project/Python') 

'''
All steps in the code is defined as a function. 
The main function is the last function and from where the input is loaded and all the other functions are called.
To read and understand the code, start from the main function, and go back to the other functions ones they are called.
'''

# Check if an input file exists
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')

# set-up the input data before training
def data_setup(gene_dataset, use_batches=False, batch_id = None):
    # Extract the cell_IDs for later use
    cells = pd.DataFrame(sc.get.obs_df(gene_dataset).index)
    cells = cells.rename(columns={0: "cell_ID"})
    # save the raw count data so that it is kept after normalization
    gene_dataset.layers['counts'] = gene_dataset.X.copy()    
    if use_batches == True:
        # define cell_IDs
        cellids = list(sc.get.obs_df(gene_dataset).index)
        # cut out the part of the cell_IDs relevand for the batch identification
        for i in range(0, len(cellids)):
            cid = []
            for z in range(0, len(batch_id)):
                cid.append(cellids[i].split('_')[batch_id[z]])
            cellids[i] = '_'.join(cid)   
        print('Use the following batches; ', sorted(set(cellids)))
        # assign numbers to each unique batch ID
        batch = list(pd.factorize(cellids)[0])        
        # add batch information to gene_dataset
        df = pd.DataFrame({'batch':batch,'cellids':cellids})        
        df.index = gene_dataset.obs.index
        gene_dataset.obs = df        
    return gene_dataset, cells

# function for cell clustering. 
# Note: This function may need to be updated for the new scVI-tools
def cluster_cells(full, outdir, cells, resolution, cellids = None):
    '''
    The "cells" here can be exchanged for full.gene_dataset.CellID to clean up the code
    '''
    latent, batch_indices, labels = full.sequential().get_latent()
    # Visualizing latent features
    post_adata = sc.AnnData(latent)
    sc.pp.neighbors(post_adata)
    sc.tl.umap(post_adata)
    sc.pl.umap(post_adata, save=True)

    #Computing clusters, resolution affects the amount of clusters
    sc.tl.leiden(post_adata, key_added="leiden_scvi", resolution=resolution)
    #Plotting clusters
    sc.pl.umap(post_adata, color=["leiden_scvi"], save=True) 
    '''
    When testing the clustering with different parameters, check the fig 
    at this stage to see if it looks reasonable. 
    Continue when all parameters are ok.
    '''
    # Move the figure to the output directory
    if os.path.exists(outdir + "/figures") == False:
        print('dir_out was created')
        os.makedirs(outdir + "/figures")
    os.rename(os.getcwd() + "/figures/umap.pdf", outdir + "/figures/umap.pdf")
    # Plot with specific parameters
    if cellids != None:
        post_adata.obs['batch'] = cellids
        sc.pl.umap(post_adata, color=["batch"], save=True)
        os.rename(os.getcwd() + "/figures/umap.pdf", outdir + "/figures/umap_batchlabels.pdf")
    if len(os.listdir('figures/')) == 0:
        os.rmdir('figures')
    #Extracting cluster labels
    clusterID = post_adata.obs["leiden_scvi"]
    cells['clusterID'] = list(clusterID)
    cells.to_csv(outdir + '/cluster_ids.csv', index=False)

# set-up data before DEG analysis with cluster information etc. 
def DEG_setup(clusts, gene_dataset, group_column, sort_column = None):
    # identify the colnames of the groupid
    groupid = clusts.columns[group_column]
    if sort_column != None:
        # identify the sortid names based on the colnames in the cluster input file
        sortid = []
        for z in range(0, len(sort_column)):
            sortid.append(clusts.columns[sort_column[z]])
        # identify the lables for sorting
        lables = np.array(clusts[[sortid[0]]].values)
        if len(sortid) > 1:
            for z in range(1, len(sortid)):
                lables = np.array(lables + '_' + clusts[[sortid[z]]].values)
        # add the groupids to the lables
        lables = np.array(lables + '_' + clusts[[groupid]].values)
        # create the list of sortids
        sortid_out = []
        for z in range(0, len(sortid)):
            sortid_out.append(np.unique(clusts[[sortid[z]]].values))
        # identify the unique lables
        lables_unique = np.unique(lables)
        # assign numerical IDs to each unique lable
        lables_num = []
        for row in range(0,len(lables)):
            lables_num.append(list(lables_unique).index(lables[row]))
        lables_num = np.array(lables_num)
        # add lable ids to the cluster dataframe       
        clusts[['lables']] = lables 
        clusts[['lables_num']] = lables_num 
        clusts.head()
        # check that the dimentions of the custer table and the expression matrix are the same,
        # and add the cluster information to previous cell information in the gene_dataset
        if (len(clusts) != len(gene_dataset.obs)):
            print('The length of the cluster information and gene_dataset is not the same.')
            exit
        gd = gene_dataset.obs
        gd['cell_ID'] = gd.index
        clusts = pd.merge(gd,clusts,left_on='cell_ID',right_on='cell_ID',how='inner')
        clusts.index = clusts['cell_ID'].values
        gene_dataset.obs = clusts
        # identify the unique groupids
        groupid = np.unique(clusts[[groupid]].values)
        return gene_dataset, sortid_out, groupid, lables, lables_unique
    elif sort_column == None:
        
        lables = np.array(clusts[[groupid]].values)
        lables_unique = np.unique(lables)
        lables_num = []
        for row in range(0,len(lables)):
            lables_num.append(list(lables_unique).index(lables[row]))
        lables_num = np.array(lables_num)
        # add lable ids to the cluster dataframe       
        clusts[['lables']] = lables 
        clusts[['lables_num']] = lables_num 
        clusts.head()

        # check that the dimentions of the custer table and the expression matrix are the same,
        # and add the cluster information to previous cell information in the gene_dataset
        if (len(clusts) != len(gene_dataset.obs)):
            print('The length of the cluster information and gene_dataset is not the same.')
            exit
        gd = gene_dataset.obs
        gd['cell_ID'] = gd.index
        clusts = pd.merge(gd,clusts,left_on='cell_ID',right_on='cell_ID',how='inner')
        clusts.index = clusts['cell_ID'].values
        gene_dataset.obs = clusts
        # identify the unique groupids
        groupid = np.unique(clusts[[groupid]].values)
        return gene_dataset, groupid, lables, lables_unique

# DEG analysis when len(groupid) >= 2
def DEG_analysis_groupid1(gene_dataset, full, couple_celltypes, lables_unique, outdir_DEG, outname):
    if len(couple_celltypes[0]) == 2:
        # create boolean arrays of which cells in the gene_dataset to calculate DEGs between
        cell_idx1 = gene_dataset.obs[['lables_num']].values == list(pd.DataFrame(couple_celltypes)[1])
        cell_idx2 = gene_dataset.obs[['lables_num']].values == list(pd.DataFrame(couple_celltypes)[0])
        # Make sure that each group contains at least three cells. If not, 
        # the resulting lists of DEGs may contain a lot of false positives. 
        if np.sum(cell_idx1) < 3 or np.sum(cell_idx2) < 3:
            print('Skip to next. Less than 3 cells in one of the groups')
        else:
            print("\nDifferential Expression A/B for cell types\nA: %s\nB: %s\n" %
                  tuple((lables_unique[pd.DataFrame(couple_celltypes)[i]] for i in [1, 0])))
            # calculate DEGs based on vanilla mode and write to out
            #'''
            #vanilla is the prefered method for calculation of DEGs.
            #This method requires FCs to be calculated separetly             
            #'''
            #de_vanilla = full.differential_expression(
            #        idx1=cell_idx1,
            #        idx2=cell_idx2,
            #        mode='vanilla'
            #)
            #print('test: deg done')
            #print('write output file: ' + outname)
            #de_vanilla.to_csv(outdir_DEG + '/vanilla_mode/' + outname)

            # calculate DEGs based on change mode and write to out
            '''
            change is an alternative method to calculate DEGs, 
            which also includes calculations of FCs. 
            '''
            de_change = full.differential_expression(
                    idx1=cell_idx1,
                    idx2=cell_idx2,
                    mode='change'
            )
            print('write output file: ' + outname)
            de_change.to_csv(outdir_DEG + '/change_mode/' + outname)                                
    elif len(couple_celltypes[0]) == 1:
        # create boolean arrays of which cells in the gene_dataset to calculate DEGs between
        cell_idx1 = gene_dataset.obs[['lables_num']].values == list(pd.DataFrame(couple_celltypes)[0])
        cell_idx2 = gene_dataset.obs[['lables_num']].values != list(pd.DataFrame(couple_celltypes)[0])
        # Make sure that each group contains at least three cells. If not, 
        # the resulting lists of DEGs may contain a lot of false positives. 
        if np.sum(cell_idx1) < 3 or np.sum(cell_idx2) < 3:
            print('Skip to next. Less than 3 cells in one of the groups')
        else:
            print("\nDifferential Expression A/B for cell types\nA: %s\nB: AllOther" %
                  tuple((lables_unique[pd.DataFrame(couple_celltypes)[i]] for i in [0])))
            # calculate DEGs based on vanilla mode and write to out
#            '''
#            vanilla is the prefered method for calculation of DEGs.
#            This method requires FCs to be calculated separetly             
#            '''
#            de_vanilla = full.differential_expression(
#                    idx1=cell_idx1,
#                    idx2=cell_idx2,
#                    mode='vanilla'
#            )
#            print('test: deg done')
#            print('write output file: ' + outname)
#            de_vanilla.to_csv(outdir_DEG + '/vanilla_mode/' + outname)

            # calculate DEGs based on change mode and write to out
            '''
            change is an alternative method to calculate DEGs, 
            which also includes calculations of FCs. 
            '''
            de_change = full.differential_expression(
                    idx1=cell_idx1,
                    idx2=cell_idx2,
                    mode='change'
            )
            print('write output file: ' + outname)
            de_change.to_csv(outdir_DEG + '/change_mode/' + outname)                                
        
    else:
        print('Skip to next. Expected two groups for comparision, found ' + str(len(couple_celltypes)))

        

def main():
    '''
    If running the script from the shell (run_scVI.sh), this part is responsible 
    for identification of the input files
    '''
    parser = argparse.ArgumentParser(description='Analyse single cell data using scVI')
    parser.add_argument('--infile', help='Relative path to the input data. The input file need to be a comma separated csv file, with cells as columns and genes as rows', required=True, type=lambda x: is_valid_file(parser, x))
    parser.add_argument('--outdir', help='The path to the cluster_analysis output directory', required=True)
    parser.add_argument('--use_batches', dest='use_batches', action='store_true', default=False, help='add if normalization over batches are needed.')
    parser.add_argument('--batch_id', nargs='+', type=int, default=0, help='Define which parts of the colnames (infile) should be used to define batches (count from 0, sep = "_"). eg cellname; "BatchID_cell1". Default = 0. If multiple sources of batch effect, input a space-separated list of integers')
    parser.add_argument('--n_epochs', type=int, default=400, help='This argument will define the amount of training and should depend on the number of cells in the input data. Recommended amount is 100 or more. Default = 400.')
    parser.add_argument('--cluster_analysis', dest='cluster_analysis', action='store_true', default=False, help='add if you want to perform the cluster analysis.')
    parser.add_argument('--resolution', type=float, default=0.5, help='This argument will define the resolution of clustering. Default = 0.5. Only used if cluster_analysis == True')
    parser.add_argument('--DEG_analysis', dest='DEG_analysis', action='store_true', default=False, help='add this command if you want to perform the differential expression analysis.')
    parser.add_argument('--clust_file', help='Relative path to the cluster data file. Required for --DEG_analysis', required=False)
    parser.add_argument('--outdir_DEG', help='The path to the DEG_analysis output directory. Required for --DEG_analysis', required=False)
    parser.add_argument('--sort_column', nargs='+', type=int, default=None, help='Column in clust_file to sort the data for DEG analysis. DEGs will be calculated for each unique variable in this column separetly. (count first column as 0). If multiple sources to sort on, input a space-separated list of integers')
    parser.add_argument('--group_column', nargs='+', type=int, default=None, help='Column in clust file containing the groups between which to identify DEGs. If the column contains more that two groups, DEGs will be calculated for each unique value vs all others. (count first column as 0)')

    args = parser.parse_args()

    '''
    If some sorting step of the cells are to be included in this script,
    check how the cell names are introduced in the results. This will probably need to
    be done in a different way so to track which cells that are removed.

    Update ideas: Is it possible to introduce the cell names in a clearer way when
    loading the input data.
    '''

    print('Initiate analysis')

    '''
    Here is the loading of the input data. This part only works if run from command line or shell script
    '''
    infile = os.getcwd() + '/' + args.infile.name
    outdir = os.getcwd() + '/' + args.outdir
    use_batches = args.use_batches
    if use_batches == True:
        batch_id = tuple(args.batch_id)
    n_epochs = args.n_epochs
    cluster_analysis = args.cluster_analysis
    DEG_analysis = args.DEG_analysis
    if cluster_analysis == True:
        resolution = args.resolution
    if DEG_analysis == True:
        clust_file = os.getcwd() + '/' + args.clust_file
        outdir_DEG = os.getcwd() + '/' + args.outdir_DEG
        if args.sort_column!=None:
            sort_column = tuple(args.sort_column)
        group_column = tuple(args.group_column)

    '''
    To test the script from the python environment (I use spyder), the data can be loaded as in the commented lines.    
    '''
#    infile = os.getcwd() + "/../data/sorted_DGEs/sorted_expression_matrix.csv"
#    outdir = os.getcwd() + "/../data/scVI_normalized/"
#    use_batches = True
#    batch_id = tuple([0,2])
#    n_epochs = 400
#    cluster_analysis = False

#    DEG_analysis = True
#    clust_file = os.getcwd() + "/../data/clusters_final_out/cluster_ids_fromOleg_06_29.csv"
#    outdir_DEG = os.getcwd() + "/../results/DEG_analysis/cluster_ids_fromOleg_06_29_CellType/marker_genes"

#    sort_column = tuple([4])
#    group_column = tuple([3])


    print('load and set-up the data')
    gene_dataset = scvi.data.read_csv(infile)
    gene_dataset = gene_dataset.transpose()

    if use_batches == True:
        gene_dataset, cells = data_setup(gene_dataset, use_batches, batch_id)
    elif use_batches == False:
        gene_dataset, cells = data_setup(gene_dataset, use_batches)

    # Normalizing the data
    sc.pp.normalize_total(gene_dataset, target_sum=1e4)
    sc.pp.log1p(gene_dataset)
    gene_dataset.raw = gene_dataset

    # set-up the anndata (scvi) variable. This is an important step for scVI training and other analyses. 
    # If changes are made to the gene_dataset, this need to be set up again. 
    if use_batches == True:
        scvi.data.setup_anndata(gene_dataset, layer="counts", batch_key="batch")
    elif use_batches == False:
        scvi.data.setup_anndata(gene_dataset, layer="counts")
    

    save_dir = os.path.join(outdir + "/full_posterior_" + str(n_epochs))
    
    # Define and train model
    '''
    This is a time consuming step, and each run will give slightly different results. 
    For that reason, the trained data is saved so that if one need to re-run some analyses, 
    it can be done without much time being consumed. 
    Remember, that the training is dependent on the input gene_dataset. If something is changed 
    in the earlier parts of this code, this training step need to be re-run. 
    '''
    if not os.path.exists(save_dir + '/model_params'):
        print('train the data')
        full = scvi.model.SCVI(gene_dataset)    
        full.train(n_epochs = n_epochs)  
        # ... after training, save your model 
        full.save(save_dir + '/model_params/')
    elif os.path.exists(save_dir + '/model_params'):
        print('load the trained data')    
        full = scvi.model.SCVI.load(save_dir + '/model_params/', gene_dataset, use_cuda=True)

    # get latent variables
    latent = full.get_latent_representation() 
    #latent_MonteCarlo_1000 = full.get_latent_representation(give_mean=False, mc_samples=1000) 
    #latent_MonteCarlo_1000 = full.get_latent_representation(mc_samples=1000) 
    gene_dataset.obsm["X_scVI"] = latent
    # get normalized matrix
    normalized = full.get_normalized_expression()
    # get latent output
    latent_out = pd.DataFrame(latent)
    latent_out.index = normalized.index
    #latent_MonteCarlo_1000_out = pd.DataFrame(latent_MonteCarlo_1000)
    #latent_MonteCarlo_1000_out.index = normalized.index
    # write to out, unless already exist
    '''
    These matrices may be needed for downstream analyses. I save them now to keep them ready.
    '''
    if not os.path.exists(outdir + '/latent_expression_matrix.csv'):
        latent_out.to_csv(outdir + '/latent_expression_matrix.csv')
 #   if not os.path.exists(outdir + '/latent_MonteCarlo_1000_expression_matrix.csv'):
 #       latent_MonteCarlo_1000_out.to_csv(outdir + '/latent_MonteCarlo_1000_expression_matrix.csv')
    if not os.path.exists(outdir + '/normalized_expression_matrix.csv'):
        normalized.to_csv(outdir + '/normalized_expression_matrix.csv')
    

    # Cluster analysis
    '''
    This part of the script can be run in order to cluster the cells. It will be nessesary 
    to have the clusters before DEG analysis, to know which groups to calculate the DEGs between. 
    Clusters can be achieved by this method, or by some other method. 
    
    Note: Thgis method is not updated yet for the new version of scVI-tools. 
    '''
    if cluster_analysis == True:
        print('Update the code for new version of scVI-tools')
        exit
        # clustering analysis
        print('Perform cluster analysis')
        if use_batches == True:
            cluster_cells(full, outdir, cells, resolution, gene_dataset.obs['cellids'])
        elif use_batches == False:
            cluster_cells(full, outdir, cells, resolution)

    # Differential expression analysis
    '''
    This part of the code is used to calculate DEGs between groups of cells, eg sick vs control. 
    '''
    if DEG_analysis == True:
        #if DEG_XvsY == True:
#        import inspect
#        test = inspect.getsource(full.differential_expression)
        if 'clust_file'  not in locals():
            raise ValueError('The clust_file is missing')
        print('set-up data for DEG analysis')
        # load the cluster information. This file will be used to define the different groups of cells. 
        clusts = pd.read_csv(clust_file, dtype=str)
        
        # Update the gene_dataset variable based on the cluster information
        if args.sort_column==None:
            gene_dataset, groupid, lables, lables_unique = DEG_setup(clusts, gene_dataset, group_column)
        else:
            gene_dataset, sortid, groupid, lables, lables_unique = DEG_setup(clusts, gene_dataset, group_column, sort_column)

        # As changes have not been made to the gene_dataset variable, set-up the anndata (scvi) variable again. 
        # Make sure that it is set-up in the same way as before training, or you will get errors while running the code. 
        if use_batches == True:
            scvi.data.setup_anndata(gene_dataset, layer="counts", batch_key="batch")
        elif use_batches == False:
            scvi.data.setup_anndata(gene_dataset, layer="counts")

        '''
        if statements are to define the arguments in the correct way dependent on how many groupids and sortids 
        you have. 
        
        The groupid is the groups between which you want to calculat4e DEGs, eg. sick and healthy (len(groupid) == 2)
        or cell type X vsw all other cell types (len(groupid) > 2)
        
        The sortid contains the information about the clusters and other information improtant to subset the data 
        before DEG analysis. Eg cell types (len(ortid) == 1), cell types and tissues (len(ortid) == 2) etc.
        
        All infromation for groupid and sortid are presented in the input cluster file, and the columns to consider 
        for these parts are included in the input variables. 
        '''
        if len(groupid) < 2:
            raise ValueError('Less than 2 groups, cannot calculate DEGs')
        if len(groupid)==2:
            # put into function
            print('Perform DEG analysis')
            if args.sort_column==None:
                '''
                No sort_column mean that DEGs will be calculated between all cells of the groupids,
                i.e. no previous subseting of clusters
                
                Note: args argument is only working when scripts is run from command line
                '''
                outname = 'DEGs_' + groupid[1] + '_vs_' + groupid[0] + '.csv'
                if os.path.exists(outdir_DEG + '/' + outname):
                    print()
                    sys.exit(outname + ' has already been produced.')
                print('find DEGs')
                print('Set populations for comparisions')
                '''
                Note: The exact pattern of groupid will be searched among the lables_unique.
                '''
                couple_celltypes = np.where(np.isin(lables_unique,groupid))                
                # DEG analysis 
                DEG_analysis_groupid1(gene_dataset, full, couple_celltypes, lables_unique, outdir_DEG, outname)

            elif len(sortid) == 1:
                '''
                For each group (eg cell type) in the sample, 
                the DEGs will be calculated between eg. sick and healthy
                '''
                for sid1 in sortid[0]:
                    sid = sid1
                    outname = 'DEGs_' + groupid[1] + '_vs_' + groupid[0] + '_' + re.sub(' ', '-', sid) + '.csv'
                    if os.path.exists(outdir_DEG + '/' + outname):
                        print(outname + ' has already been produced. Skip to next.')
                        continue
                    print('find DEGs for ' + sid)
                    print('Set populations for comparisions')
                    '''
                    Note: The exact pattern formed by sid + '_' + groupid will be searched among the lables_unique.
                    Therefore, make sure the input files looks as the example data, or this part may need to be updated
                    '''
                    couple_celltypes = np.where(np.isin(lables_unique,sid + '_' + groupid))
                    # DEG analysis                    
                    DEG_analysis_groupid1(gene_dataset, full, couple_celltypes, lables_unique, outdir_DEG, outname)

            elif len(sortid) == 2:
                '''
                For each group combination (eg cell type X in tissue Y) in the sample, 
                the DEGs will be calculated between eg. sick and healthy
                '''
                for sid1 in sortid[0]:
#                    sid1 = sortid[0][8]
                    for sid2 in sortid[1]:
#                        sid2 = sortid[1][0]
                        sid = sid1 + '_' + sid2
                        outname = 'DEGs_' + groupid[1] + '_vs_' + groupid[0] + '_' + re.sub(' ', '-', sid) + '.csv'
                        if os.path.exists(outdir_DEG + '/' + outname):
                            print(outname + ' has already been produced. Skip to next.')
                            continue
                        print('find DEGs for ' + sid)
                        print('Set populations for comparisions')
                        '''
                        Note: The exact pattern formed by sid + '_' + groupid will be searched among the lables_unique.
                        Therefore, make sure the input files looks as the example data, or this part may need to be updated
                        '''
                        couple_celltypes = np.where(np.isin(lables_unique,sid + '_' + groupid))                        
                        # DEG analysis 
                        DEG_analysis_groupid1(gene_dataset, full, couple_celltypes, lables_unique, outdir_DEG, outname)
                        
            elif len(sortid) > 2:
                raise ValueError('Sort by more than 2 categories, code not developed for this yet')

        elif len(groupid) > 2:
            '''
            The script will compare group(i) vs group(allOther)            
            '''
            print('compute DEGs between celltype X and all other cell types for: ' + str(groupid))

            if args.sort_column==None:
                for group in groupid:
#                    group = groupid[0]
                    outname = 'DEGs_' + re.sub(' ', '-', group) + '_vs_AllOthers.csv'
                    if os.path.exists(outdir_DEG + '/' + outname):
                        print(outname + ' has already been produced. Skip to next.')
                        continue
                    print('find DEGs for ' + group + ' vs all others')
                    print('Set populations for comparisions')
                    '''
                    Note: The exact pattern formed by sid + '_' + groupid will be searched among the lables_unique.
                    Therefore, make sure the input files looks as the example data, or this part may need to be updated
                    '''
                    couple_celltypes = np.where(np.isin(lables_unique,group))                        
                    # DEG analysis 
                    DEG_analysis_groupid1(gene_dataset, full, couple_celltypes, lables_unique, outdir_DEG, outname)



            elif len(sortid) == 2:
                print('Update the code for new version of scVI-tools')
                exit
                for sid1 in sortid[0]:
#                    sid1 = sortid[0][0]
                    for sid2 in sortid[1]:
#                        sid2 = sortid[1][0]
                        sid = sid1 + '_' + sid2
                        print('find DEGs for ' + sid)
                        print('Set populations for comparisions')
                        #couple_celltypes = np.where(np.isin(lables_unique,sid + '_' + groupid))
                        #couple_celltypes_names = lables_unique[couple_celltypes]

                        # repeate the following line for len(sortid) amount of times
                        couple_celltypes_names = lables_unique[np.where([sid1 in i for i in lables_unique])]                        
                        couple_celltypes_names = couple_celltypes_names[np.where([sid2 in i for i in couple_celltypes_names])]
                        couple_celltypes = np.where(np.isin(lables_unique,couple_celltypes_names))                        
                        
                        for i in range(0, len(couple_celltypes_names)):
#                            i = 0
                            couple_celltypes_names
                            outname = 'DEGs_' + re.sub(' ', '-', couple_celltypes_names[i]) + '_vs_allOtherCells.csv'
                            if os.path.exists(outdir_DEG + '/' + outname):
                                print(outname + ' has already been produced. Skip to next.')
                                continue
                            cell_idx1 = gene_dataset.labels.ravel() == list(pd.DataFrame(couple_celltypes))[i]
                            couple_celltypes_list = list(pd.DataFrame(couple_celltypes))
                            couple_celltypes_list.remove(i)
                            cell_idx2 = np.isin(gene_dataset.labels.ravel(), couple_celltypes_list)
                            if np.sum(cell_idx1) < 3 or np.sum(cell_idx2) < 3:
                                print('Skip to next. Less than 3 cells in one of the groups')
                                continue
                            print("\nDifferential Expression for cell type vs all other cell types\n%s\n" %
                                 tuple((lables_unique[pd.DataFrame(couple_celltypes)[i]])))
#                            de_vanilla = full.differential_expression_score(
#                                cell_idx1,
#                                cell_idx2,
#                                use_observed_batches=True,
#                            )
#                            print('write output file: ' + outname)
#                            de_vanilla.to_csv(outdir_DEG + '/' + outname)
#
                            de_change = full.differential_expression_score(
                                cell_idx1,
                                cell_idx2,
                                use_observed_batches=True,
                                mode='change',
                            )
                            print('write output file: ' + outname)
                            de_change.to_csv(outdir_DEG + '/FCs/' + outname)



if __name__ == '__main__':
   main()

