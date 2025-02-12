"""
This Module is to run the modifinder algorithm on the give input data.
"""

from modifinder import ModiFinder, Compound, BasicEvaluationEngine
import pandas as pd
import pickle

def load_Compound_from_cache(data, cache = None, **kwargs):
    """
    Load the given data from the cache.
    """
    if isinstance(data, str):
        if cache is not None and data in cache:
            return Compound(cache[data], id = data, **kwargs)
    
    return Compound(data, **kwargs)

def run_single(match_index, network = None, networkUnknowns = None, unknown_compound = None, known_compounds = None, 
               helpers = None, cached_compounds = None, **kwargs):
    """
    Run the modifinder algorithm on the given input data and return the result as a dictionary.
    """
    
    evaluation_engine = BasicEvaluationEngine(**kwargs)
    final_result = [dict() for _ in range(len(known_compounds))]
    try:
        if network is not None:
            mf = ModiFinder(network = network, networkUnknowns = networkUnknowns, **kwargs)
        else:
            mf = ModiFinder()
            unknown_compound = load_Compound_from_cache(unknown_compound, cached_compounds, **kwargs)
            unknown_id = mf.add_unknown(unknown_compound, **kwargs)
            for known_index, known_compound in enumerate(known_compounds):
                known_compound = load_Compound_from_cache(known_compound, cached_compounds, **kwargs)
                mf.add_neighbor(known_compound, unknown_id)
                for helper in helpers[known_index]:
                    helper = load_Compound_from_cache(helper, cached_compounds, **kwargs)
                    mf.add_neighbor(helper, known_compound.id)
                    
        mf.re_annotate(mf.annotationEngine, **kwargs)
        mf.re_align(mf.alignmentEngine, **kwargs)
        
        knowns = mf.get_neighbors(node_id=unknown_id, **kwargs)
        
    except Exception as err:
        final_result = [{"match_index": match_index, "error": str(err), match_index: match_index} for _ in range(len(known_compounds))]
        try:
            for known_index, node in enumerate(known_compounds):
                if isinstance(node, str):
                    final_result[known_index]["known_id"] = node
                if isinstance(unknown_compound, str):
                    final_result[known_index]["unknown_id"] = unknown_compound
                elif unknown_compound is not None:
                    final_result[known_index]["unknown_id"] = unknown_compound.id
        except:
            pass
        # raise err
        return final_result
    
    for known_index, node in enumerate(knowns):
        try:
            probs = mf.generate_probabilities(unknown_id=unknown_id, known_id=node, **kwargs)
            data = mf.get_meta_data(unknown_id=unknown_id, known_id=node)
            result = {
                "match_index": match_index,
                "unknown_id": unknown_id,
                "known_id": node,
            }
            result.update(data)
            
            if unknown_compound.structure:
                for method in ["is_max", "average_distance"]:
                    known_compound = mf.network.nodes[node]["compound"]
                    evaluation_result = evaluation_engine.evaluate_single(known_compound.structure, unknown_compound.structure,
                                                                            probs, evaluation_method=method, **kwargs)
                    result[method] = evaluation_result
        
            final_result[known_index] = result
                    
        except Exception as err:
            final_result[known_index] = {"error": str(err), "match_index": match_index, "unknown_id": unknown_id, "known_id": node}
            # raise err
    return final_result
    

def run_batch(matches, output_dir, save_pickle = True, save_csv = True, cached_data = None) -> pd.DataFrame:
    """
    run the modifinder on the batch of matches
    """
    result = []
    for match in matches:
        result.extend(run_single(**match, cached_compounds = cached_data))
    
    df = pd.DataFrame(result)
    
    if save_pickle:
        with open(output_dir + ".pkl", "wb") as f:
            pickle.dump(df, f)
    
    if save_csv:
        df.to_csv(output_dir + ".csv", index = False)
    
    
    

    