
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def resample_med(df1, df2, col, essential_genes, direction, n=10000, graph=False, noisy=True,figName='resampling_distribution.png'):

    # inputs: df1 = DataFrame of candidate genes, needs column of gene names named 'gene'
    #         df2 = DataFrame of all genes to resample from, needs column of gene names named 'gene'
    #         col = name of df column with statistic to resample
    #         essential = list of essential genes
    #         direction = ['less_than', 'greater']: test whether statistic is less than or greater than random samples
    #         n = number of random samples to generate
    #         graph = show KDE plot of resampling distribution
    #         noisy = print statistics

    
    # output: p value (proportion of random samples with statistic ['less_than', 'greater_than'] or equal to candidates


    actual_med = df1[col].median()
    
    essential_count = len([i for i in df1['gene'] if i in essential_genes])
    nonessential_count = len(df1[col]) - essential_count

    essential_df = df2[df2['gene'].isin(essential_genes)]
    nonessential_df =  df2[~df2['gene'].isin(essential_genes)]

    sample_meds = []
    results = 0

    for _ in range(n):
        s = random_sample_med(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_med, col)
        results += s[0]
        sample_meds += [s[1]]
    
    if graph:  
        #sns.displot(sample_meds, kind='kde') #new seaborn
        sns.distplot(sample_meds, kde=True) #old seaborn
        plt.axvline(x=actual_med, color='r', label = 'True median')
        plt.xlabel('Sample medians')
        plt.ylabel('Frequency')
        plt.title(col + ' resampling distribution, p= ' + str(results/n))
        plt.savefig(figName)

    if noisy:
        print('candidate gene median {}: {}'.format(col, actual_med))
        print('essential count: {}; nonessential_count: {}'.format(essential_count, nonessential_count))
        print('resampling pool size: {}'.format(df2.shape[0]))
        print('p = {}'.format(results/n))
    
    return results/n       


def random_sample_med(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_med, col):

        sample = essential_df.sample(n=essential_count)
        sample = sample.append(nonessential_df.sample(n=nonessential_count))
        sample_med = sample[col].median()

        rv = 0
        if direction == 'less_than':
            if sample_med <= actual_med:
                rv += 1

        elif direction == 'greater_than':
            if sample_med >= actual_med:
                rv += 1

        return rv, sample_med


### as above but for means
### TO DO: combine resample_mean and resample_median functions into one


def resample_mean(df1, df2, col, essential_genes, direction, n=10000, graph=False, noisy=True):

    # inputs: df1 = DataFrame of candidate genes, needs column of gene names namean 'gene'
    #         df2 = DataFrame of all genes to resample from, needs column of gene names namean 'gene'
    #         col = name of df column with statistic to resample
    #         essential = list of essential genes
    #         direction = ['less_than', 'greater']: test whether statistic is less than or greater than random samples
    #         n = number of random samples to generate
    #         graph = show KDE plot of resampling distribution
    #         noisy = print statistics

    
    # output: p value (proportion of random samples with statistic ['less_than', 'greater_than'] or equal to candidates


    actual_mean = df1[col].mean()
    
    essential_count = len([i for i in df1['gene'] if i in essential_genes])
    nonessential_count = len(df1[col]) - essential_count
    
    essential_df = df2[df2['gene'].isin(essential_genes)]
    nonessential_df =  df2[~df2['gene'].isin(essential_genes)]

    sample_means = []
    results = 0

    for _ in range(n):
        s = random_sample_mean(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_mean, col)
        results += s[0]
        sample_means += [s[1]]
    
    if graph:  
        sns.displot(sample_means, kind='kde')
        plt.axvline(x=actual_mean, color='r', label = 'True mean')
        plt.xlabel('Sample means')
        plt.ylabel('Frequency')
        plt.title(col + ' resampling distribution, p= ' + str(results/n))

    if noisy:
        print('candidate gene mean {}: {}'.format(col, actual_mean))
        print('essential count: {}; nonessential_count: {}'.format(essential_count, nonessential_count))
        print('resampling pool size: {}'.format(df2.shape[0]))
        print('p = {}'.format(results/n))
    
    return results/n       


def random_sample_mean(essential_df, essential_count, nonessential_df, nonessential_count, direction, actual_mean, col):

        sample = essential_df.sample(n=essential_count)
        sample = sample.append(nonessential_df.sample(n=nonessential_count))
        sample_mean = sample[col].mean()

        rv = 0
        if direction == 'less_than':
            if sample_mean <= actual_mean:
                rv += 1

        elif direction == 'greater_than':
            if sample_mean >= actual_mean:
                rv += 1

        return rv, sample_mean
