from modifinder.engines.Abtracts import PredictionEngine
from depricated.Compound import Compound
import numpy as np

class BasicPredictionEngine(PredictionEngine):
    def __init__(self, **kwargs):
        """
        Initializes the basic prediction engine

        Parameters:
            :kwargs: additional arguments
        """
        self.args = kwargs
        pass

    def predict(self, network, **kwargs):
        """
        Predicts the modifications in the network

        Parameters:
            :network (Network): the network to be predicted
            :kwargs: additional arguments
        """
        pass

    def predict_single(self, compound, modify_compound=True, **kwargs):
        """
        Predicts the modifications in a single compound

        Parameters:
            :compound (Compound): the compound to be predicted
            :modify_compound (bool): whether to modify the passed compound with the predictions
            :kwargs: additional arguments
        
        Returns:
            :List[Tuple[int, List[str]]]: a list of tuples with the indices of the peaks as keys and the list of predictions as values
        """
        pass

def find_existance(compound: Compound, peakids, fragment_engine_instance):
        """
        For the given peaks, finds the fragments that each atom contributes to
        input:
        peakids: list of peak ids
        output:
        existance: list of dicts for each atom, each dict contains the peak ids as keys and the fragments as values
        """
        existance = [dict() for i in range(len(compound.structure.GetAtoms()))]
        for peak in peakids:
            for fragment in compound.peak_fragments_map[peak]:
                hitAtoms = fragment_engine_instance.get_fragment_info(fragment, 0)[1]
                for atom in hitAtoms:
                    if peak not in existance[atom]:
                        existance[atom][peak] = []
                    existance[atom][peak].append(fragment)
        return existance
    
def calculate_contribution_atom_in_peak(compound, atom, peak, existance_data, CI = False, CPA = True, CFA = True):
    contribution = 0
    if peak not in existance_data[atom]:
        return contribution
    
    intensity_factor = 1
    atom_peak_ambiguity_factor = 1
    fragment_ambiguity_factor = 1

    if CI:
        intensity_factor = compound.peaks[peak][1]
    if CPA:
        atom_peak_ambiguity_factor = 1/len(compound.peak_fragments_map[peak])

    for frag in existance_data[atom][peak]:
        if CFA:
            # fragment_ambiguity_factor = 1 - len(compound.fragments.get_fragment_info(frag, 0)[1])/len(compound.structure.GetAtoms())
            fragment_ambiguity_factor = 1/len(compound.fragments.get_fragment_info(frag, 0)[1])
        
        contribution += intensity_factor * atom_peak_ambiguity_factor * fragment_ambiguity_factor
    
    return contribution

def calculate_contributions(compound, peakids, fragment_engine_instance, CI = False, CPA = True, CFA = True, CPE = True):
    """ 
    input:
    peakids: list of peak ids
    CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
    CPA: (Consider_Peak_Ambiguity) bool, if True, the peak ambiguity (number of fragments assigned to a peak) is considered (default: True)
    CFA: (Consider_Fragment_Ambiguity) bool, if True, the fragment ambiguity (number of atoms in fragment) is considered (default: True)
    CPA: (Consider_Peak_Entropy) bool, if True, the peak entropy (how ambiguis the fragments are) is considered (default: True
    """
    num_atoms = len(compound.structure.GetAtoms())
    existance_data = find_existance(compound, peakids, fragment_engine_instance)
    contributions = [0 for i in range(num_atoms)]
    peak_atom_contributions = np.zeros((len(peakids), num_atoms))
    for i, peak in enumerate(peakids):
        for atom in range(num_atoms):
            peak_atom_contributions[i][atom] = calculate_contribution_atom_in_peak(compound, atom, peak, existance_data, CI=CI, CPA=CPA, CFA=CFA)
    
    if CPE:
        peak_entropies = np.zeros(len(peakids))
        for i in range(len(peakids)):
            peak_entropies[i] = 1 - utils.entropy(peak_atom_contributions[i])
        
        # peak_entropies = peak_entropies / np.max(peak_entropies)
    else:
        peak_entropies = np.ones(len(peakids))
        
    
    for i in range(num_atoms):
        for j in range(len(peakids)):
            contributions[i] += peak_atom_contributions[j][i] * peak_entropies[j]
    
    
    # for atom in range(len(existance_data)):
    #     for peak in existance_data[atom]:
    #         contributions[atom] += self.calculate_contribution_atom_in_peak(atom, peak, existance_data, CI=CI, CPA=CPA, CFA=CFA)
    
    return contributions


def generate_probabilities(self, shifted_only = False, CI = False, CPA = True, CFA = True, CPE = True, method = "old"):
    """"
    Generate the probabilities for each atom to be the modification site.
    input:
    shifted_only: bool, if True, only the shifted peaks are considered
    CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
    CPA: (Consider_Peak_Ambiguity) bool, if True, the peak ambiguity (number of fragments assigned to a peak) is considered when calculating the contribution (default: True)
    CFA: (Consider_Fragment_Ambiguity) bool, if True, the fragment ambiguity (number of atoms in fragment) is considered (default: True)
    CPE: (Consider_Peak_Entropy) bool, if True, the peak entropy (how ambiguis the fragments are) is considered (default: True)
    method: str, the method to generate the probabilities (default: "old")
    """
    
    if method == "random_choice":
        n = len(self.main_compound.structure.GetAtoms())
        probabilities = np.zeros(n)
        random_choice = np.random.choice(n)
        probabilities[random_choice] = 1
        return probabilities
    elif method == "random_distribution":
        probabilities = np.random.rand(len(self.main_compound.structure.GetAtoms()))
        probabilities = probabilities / np.sum(probabilities)
        return probabilities
    elif method == "all_equal":
        probabilities = np.ones(len(self.main_compound.structure.GetAtoms()))
        probabilities = probabilities / np.sum(probabilities)
        return probabilities
    elif method == "random_skewed":
        target =  np.random.choice(len(self.main_compound.structure.GetAtoms()))
        probabilities = np.random.rand(len(self.main_compound.structure.GetAtoms()))
        probabilities = probabilities * self.main_compound.distances[target]
        probabilities = probabilities - np.min(probabilities)
        probabilities = probabilities / np.sum(probabilities)
        return probabilities
    elif method.startswith("multiple"):
        # set seed to be 0
        np.random.seed(0)
        res = []
        method_name = method[9:]
        for i in range(50):
            res.append(self.generate_probabilities(method=method_name))
        return res
    elif method == "multiply":
        probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
        for atom in range(len(self.main_compound.structure.GetAtoms())):
            positive_contributions = 1
            for peak in self.shifted:
                count = 0
                for fragment in self.main_compound.peak_fragments_map[peak[0]]:
                    if atom in self.main_compound.fragments.get_fragment_info(fragment, 0)[1]:
                        count += 1
                eps = 0
                if CPA:
                    eps = 1
                
                if (len(self.main_compound.peak_fragments_map[peak[0]]) + eps) == 0:
                    positive_contributions *= 0
                else:
                    positive_contributions *= (count + eps)/(len(self.main_compound.peak_fragments_map[peak[0]]) + eps)
            

            if shifted_only:
                negative_contributions = 0
            else:
                # TODO: implement this
                negative_contributions = 1
            
            probabilities[atom] = positive_contributions - negative_contributions
        
    else:
        s_peakids = self.get_main_shifted_index()
        positive_contributions = self.calculate_contributions(s_peakids, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        if not shifted_only:
            u_peakids = self.get_main_unshifted_index()
            negative_contributions = self.calculate_contributions(u_peakids, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        else:
            negative_contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        
        probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
        for i in range(len(positive_contributions)):
            probabilities[i] = positive_contributions[i] - negative_contributions[i]
        
        if np.min(probabilities) < 0:
            probabilities = probabilities - np.min(probabilities)
            # mask only the atoms with positive contribution
            # probabilities = [probabilities[i] if positive_contributions[i] > 0 else 0 for i in range(len(probabilities))]
        if np.sum(probabilities) > 0:
            probabilities = probabilities / np.sum(probabilities)

        probabilities = Calc_Scores.power_prob(probabilities)

    return probabilities