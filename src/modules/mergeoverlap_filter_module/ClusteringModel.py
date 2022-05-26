"""
Code for unsupervised learning for classifying false positive
classifications is below. 

Classes
    1. TrueGenomeFinder - abstract class for defining methods for inheriting classes
                        (template pattern)
    2. ClusteringModel - this model is used to classify a genome vetor as either being truly in
                         the sample or not.

Methods
    1. get_true_positive - used to find which cluster contains true positives
    2. get_filtered_genomes - Given a genome truly in the sample, this method
                              finds all other genomes within the same cluster
"""
from sklearn.datasets import make_classification
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from abc import ABC, abstractmethod
from typing import Dict, Union, List, Tuple, Optional


class TrueGenomeFinder(ABC):
    """ This is the abstract class for the unsupervised genome finder """

    @property
    @abstractmethod
    def model_name(self):
        """ getter for model name (encapsulation) """
        pass
    
    @property
    @abstractmethod
    def model(self):
        """ This is the model for the unsupervised classification """
        pass
    
    @abstractmethod
    def model_predict(self, x_vector):
        """
        Sine the model clusters data, there's not necessarily 
        a training set needed - it's practically just binary
        clustering. 
        """
        pass

class KMeansClustering(TrueGenomeFinder):
    """ This is the Kmeans clustering model class for the unsupervised genome finder """
    
    def __init__(self):
        self._model = KMeans(n_clusters=2)
        self._model_name = "KMeans Clustering"

    @property
    def model(self):
        """ This is the model for the unsupervised classification """
        return self._model

    @property
    def model_name(self):
        """ This is the model for the unsupervised classification """
        return self._model_name
    
    def model_predict(self, x_vector):
        """
        Sine the model clusters data, there's not necessarily 
        a training set needed - it's practically just binary
        clustering. 
        """
        FIT_MODEL = self.model.fit(x_vector)
        return FIT_MODEL.predict(x_vector)

class GMM(TrueGenomeFinder):
    """ This is the GMM clustering model class for the unsupervised genome finder """
    
    def __init__(self):
        self._model = GaussianMixture(n_components=2)
        self._model_name = "Gassiaun Mixture Model"

    @property
    def model(self):
        """ This is the model for the unsupervised classification """
        return self._model

    @property
    def model_name(self):
        """ This is the model for the unsupervised classification """
        return self._model_name
    
    def model_predict(self, x_vector):
        """
        Sine the model clusters data, there's not necessarily 
        a training set needed - it's practically just binary
        clustering. 
        """
        FIT_MODEL = self.model.fit(x_vector)
        return FIT_MODEL.predict(x_vector)

def get_true_positive(name_list, x_vector) -> Optional[str]:
    """
    Description:
        The goal for this method is to find a genome in the list of 
        genomes that is likely to be a part of the true cluster of genomes
        in the sample. An assumption used here is that the best genome will
        have a high mapq average and a high mergeoverlap percentage. 
    Input:
        1. name_list
        2. x_vector
                [mapq, overlap percentage]
    Output:
        1. High probability candidate name (str)
    """
    max_score = 0
    high_prob_genome = None
    for row_number, values in enumerate(x_vector):
        mapq = values[0]
        overlap_percentage = values[1]
        score = (mapq / 60) + overlap_percentage
        print(score)
        if score > max_score:
            max_score = score
            high_prob_genome = name_list[row_number]
    return high_prob_genome
    

def get_filtered_genomes(true_positive_name, model,  name_list):
    """
    Description:
        This method takes in a genome with the highest confidence 
        of being true, and collects all other genomes in the same
        bin. 
    Input:
        1. true_positive_name (str) - a name of a genome known to be true.
        2. model - a vector of 0s and 1s inidicating the binary clusters
                    NOTE: should have same order as the name list.
        3. name_list (List[str]) - 
    """
    cluster = None
    for index, name in enumerate(name_list):
        if name == true_positive_name: 
            cluster = model[index]
    return [name_list[index] for index, val in enumerate(model) if val == cluster]

