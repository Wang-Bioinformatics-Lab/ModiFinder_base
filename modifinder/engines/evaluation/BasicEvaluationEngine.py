from modifinder.engines.Abtracts import EvaluationEngine
from rdkit import Chem
from modifinder import exceptions
import modifinder.utilities.mol_utils as mf_mu
import numpy as np

def is_max_neighbor(G, probabilities, true_index, neighbor_count):
    max_indices = np.where(probabilities == probabilities.max())[0]
    for i in max_indices:
        if G[true_index, i] <= neighbor_count:
            return 1
    return 0

def is_max(G, probabilities, true_index):
    if probabilities[true_index] == max(probabilities):
        # find how many times the max value appears
        count = 0
        for i in range(len(probabilities)):
            if probabilities[i] == max(probabilities):
                count += 1

        return 1/count
    else:
        return 0
    
def dist_from_max(G, probabilities, true_index):
    min_dist = 100000
    max_val = max(probabilities)
    graph_diameter = np.amax(G)
    for i in range(len(probabilities)):
        if probabilities[i] == max_val:
            min_dist = min(min_dist, G[true_index, i])
    return float(min_dist/(graph_diameter))

def proximity(G, probabilities, true_index):
    min_dist = 100000
    max_val = max(probabilities)
    graph_diameter = np.amax(G)
    for i in range(len(probabilities)):
        if probabilities[i] == max_val:
            min_dist = min(min_dist, G[true_index, i])
    return (graph_diameter - min_dist)/graph_diameter

def average_dist_from_max(G, probabilities, true_index):
    max_val = max(probabilities)
    graph_diameter = np.amax(G)
    dists = 0
    count = 0
    for i in range(len(probabilities)):
        if probabilities[i] == max_val:
            value = G[true_index,i]/(graph_diameter)
            dists += value * probabilities[i]
            count += probabilities[i]
    if count == 0:
        return 0
    return float(dists/count)

def average_distance(G, probabilities, true_index):
    eps = 0.000001
    graph_diameter = np.amax(G)
    dists = 0
    count = 0
    for i in range(len(probabilities)):
        value = np.exp(-G[true_index, i]/(graph_diameter-G[true_index, i]+eps))
        # value  = (graph_diameter - G[true_index, i])/graph_diameter
        dists += value * probabilities[i]
        count += probabilities[i]
    if count == 0:
        return 0
    return float(dists/count)

def ranking_loss(G, probabilities, modificationSiteIdx):
    # find how far the index of the true modification site is from the max probability
    # sort the probabilities and keep the indices
    sorted_indices = np.argsort(probabilities)

    # reverse the indices
    sorted_indices = sorted_indices[::-1]

    # find the index of the true modification site
    true_index = np.where(sorted_indices == modificationSiteIdx)[0][0]

    # return the ranking loss
    return 1 - true_index/len(probabilities)


def sorted_rank(G, probabilities, modificationSiteIdx):
    # find how far the index of the true modification site is from the max probability
    # sort the probabilities and keep the indices
    sorted_indices = np.argsort(probabilities)

    # find the index of the true modification site
    true_index = np.where(sorted_indices == modificationSiteIdx)[0][0]

    # return the ranking loss
    return true_index/(len(probabilities)-1)


class BasicEvaluationEngine(EvaluationEngine):
    """Class to evaluate the probabilities based on the basic evaluation engines explained in the paper
    
        Parameters
        ----------
            default_method : str or None (default = "average_distance")
                The default evaluation method to use. If None, the user must specify the method.
                List of methods:
                
                    - **is_max** : returns 1 if the true modification site is the atom with the highest probability, 0 otherwise

                    - **proximity** : returns the proximity to the atom with the highest probability

                    - **average_distance** : returns the average distance to the atom with the highest probability normalized by the graph diameter

                    - **sorted_rank** : returns the rank of the true modification site in the sorted probabilities
                    
                    - **is_max_neighbor** : returns 1 if the true modification site is a neighbor of the atom with the highest probability, 0 otherwise

            kwargs : dict
                Additional arguments.
    """
    def __init__(self, default_method = "average_distance", **kwargs):
        """Initializes the BasicEvaluationEngine

        """
        self.default_method = default_method
        super().__init__(**kwargs)
    
    # TODO: Implement the evaluate method
    def evaluate(self, network, **kwargs):
        """Evaluates the network.

            .. warning::
                Not implemented yet.

            Parameters
            ----------
            network : type
                Description of the `network` parameter.
            **kwargs
                Additional keyword arguments.
        """
        raise exceptions.ModiFinderNotImplementedError("Method not implemented")
    
    def evaluate_single(self, known_compound_structure, unknonwn_compound_structure, probabilities, evaluation_method = None, **kwargs) -> float:
        """ Evaluates the probabilities based on the evaluation method
        
        Parameters
        ----------
            known_compound_structure : rdkit Chem.Mol
                The known compound structure
            
            unknonwn_compound_structure : rdkit Chem.Mol
                The unknown compound structure
            
            probabilities : List[float]
                The probabilities of each atom being the modification site
            
            evaluation_method : str
                The evaluation method to use. If None, the default method will be used.
            
            kwargs : dict
                additional arguments
        
        Returns
        -------
            float
                The evaluation score based on the evaluation method. The score is between 0 and 1.
            
        """
        
        method = evaluation_method if evaluation_method else self.default_method
        G = Chem.rdmolops.GetDistanceMatrix(known_compound_structure)
        
        true_modification_site = mf_mu.get_modification_nodes(known_compound_structure, unknonwn_compound_structure)
        if len(true_modification_site) != 1:
            raise exceptions.ModiFinderError("Only one modification site is allowed for this evaluation method")
        
        true_index = true_modification_site[0]
        
        if method == "is_max":
            return is_max(G, probabilities, true_index)
        elif method == "dist_from_max":
            return dist_from_max(G, probabilities, true_index)
        elif method == "proximity":
            return proximity(G, probabilities, true_index)
        elif method == "average_dist_from_max":
            return average_dist_from_max(G, probabilities, true_index)
        elif method == "average_distance":
            return average_distance(G, probabilities, true_index)
        elif method == "ranking_loss":
            return ranking_loss(G, probabilities, true_index)
        elif method == "sorted_rank":
            return sorted_rank(G, probabilities, true_index)
        elif method == "is_max_neighbor":
            if "neighbor_count" not in kwargs:
                raise exceptions.ModiFinderError("The neighbor count must be specified for this evaluation method")
            return is_max_neighbor(G, probabilities, true_index, kwargs["neighbor_count"])
        else:
            raise exceptions.ModiFinderError("Evaluation method not recognized")
        
        
        
        
        