"""
This Module is to run the modifinder algorithm on the give input data.
"""

from modifinder import ModiFinder, Compound
from modifinder.engines.evaluation.BasicEvaluationEngine import BasicEvaluationEngine
from modifinder.utilities import visualizer as mf_vis
from modifinder.utilities.general_utils import entropy
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import concurrent.futures
from tqdm import tqdm
import os
import argparse

def chunkify(lst, chunk_size):
    """Yield successive chunk_size chunks from lst."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def load_Compound_from_cache(data, cache = None, **kwargs):
    """
    Load the given data from the cache.
    """
    if isinstance(data, str):
        if cache is not None and data in cache:
            return Compound(cache[data], id = data, **kwargs)
    
    return Compound(data, **kwargs)

def run_single(match_index, network = None, networkUnknowns = None, unknown_compound = None, known_compounds = None, 
               helpers = None, cached_compounds = None, match_meta = None, output_dir = None, images_name = None, **kwargs):
    """
    Run the modifinder algorithm on the given input data and return the result as a dictionary.
    """
    evaluation_engine = BasicEvaluationEngine(**kwargs)
    if known_compounds is None:
        if match_meta is not None:
            for match in match_meta:
                match["match_index"] = match_index
            return match_meta
        else:
            raise ValueError("No known compounds are provided.")
    final_result = [dict() for _ in range(len(known_compounds))]
    try:
        if network is not None:
            mf = ModiFinder(network = network, networkUnknowns = networkUnknowns, **kwargs)
        else:
            mf = ModiFinder(**kwargs)
            unknown_compound = load_Compound_from_cache(unknown_compound, cached_compounds, **kwargs)
            if unknown_compound is None:
                raise ValueError("Unknown compound is not found.")
            unknown_id = mf.add_unknown(unknown_compound, **kwargs)
            for known_index, known_compound in enumerate(known_compounds):
                known_compound = load_Compound_from_cache(known_compound, cached_compounds, **kwargs)
                if known_compound is None:
                    raise ValueError(f"Known compound {known_index} is not found.")
                mf.add_neighbor(known_compound, unknown_id)
                if helpers is not None and known_index in helpers:
                    for helper in helpers[known_index]:
                        helper = load_Compound_from_cache(helper, cached_compounds, **kwargs)
                        mf.add_neighbor(helper, known_compound.id)
                    
        mf.re_annotate(mf.annotationEngine, **kwargs)
        mf.re_align(mf.alignmentEngine, **kwargs)
        
        knowns = mf.get_neighbors(node_id=unknown_id)
        
    except Exception as err:
        final_result = [{"match_index": match_index, "error": "error in creating the network. "  + str(err)} for _ in range(len(known_compounds))]
        for known_index, node in enumerate(known_compounds):
            try:
                try:
                    if isinstance(node, str):
                        final_result[known_index]["known_id"] = node
                    elif node is not None:
                        final_result[known_index]["known_id"] = node.id
                except Exception:
                    final_result[known_index]['known_id'] = str(node)
                    pass
                try:
                    if isinstance(unknown_compound, str):
                        final_result[known_index]["unknown_id"] = unknown_compound
                    elif unknown_compound is not None:
                        final_result[known_index]["unknown_id"] = unknown_compound.id
                except Exception:
                    final_result[known_index]['unknown_id'] = str(unknown_compound)
                    pass
            except Exception:
                pass
            
            try:
                final_result[known_index].update(match_meta[known_index])
            except Exception:
                pass
        return final_result
    
    for known_index, node in enumerate(knowns):
        try:
            result = {
                "match_index": match_index,
                "unknown_id": unknown_id,
                "known_id": node,
            }
            known_compound = mf.network.nodes[node]["compound"]
            if abs(known_compound.spectrum.precursor_mz - unknown_compound.spectrum.precursor_mz) < 0.1:
                # exact match
                probs = [0] * len(known_compounds)
                result["entropy"] = None
            else:
                probs = mf.generate_probabilities(unknown_id=unknown_id, known_id=node, **kwargs)
                result["entropy"] = entropy(probs)
            data = mf.get_meta_data(unknown_id=unknown_id, known_id=node)
            try:
                if images_name is not None:
                    if isinstance(images_name, str):
                        data['image_path'] = os.path.join(output_dir, images_name + str(match_index) + ".png")
                    else:
                        data['image_path'] = os.path.join(output_dir, images_name[known_index])
                    png = mf_vis.draw_molecule_heatmap(known_compound.structure, probs, show_labels=True, shrink_labels=True, show_legend=False)
                    plt.imsave(data['image_path'], png)
            except Exception as err:
                data['image_path'] = str(err)
                pass
            result.update(data)
            if unknown_compound.structure:
                for method in ["is_max", "average_distance"]:
                    evaluation_result = evaluation_engine.evaluate_single(known_compound.structure, unknown_compound.structure,
                                                                            probs, evaluation_method=method, **kwargs)
                    result[method] = evaluation_result
            
            # add entropy of the probabilities
            
            result["error"] = "No Issues"
            final_result[known_index] = result
                    
        except Exception as err:
            final_result[known_index] = {"error": "error when making prediction. " + str(err), "match_index": match_index, "unknown_id": unknown_id, "known_id": node}
            # raise err
    
    try:
        for i in range(len(known_compounds)):
            try:
                final_result[i].update(match_meta[i])
            except Exception:
                pass
    except Exception:
        pass
    return final_result
    

def run_batch(matches, output_dir, file_name, save_pickle = True, save_csv = True, cached_data = None) -> list:
    """
    run the modifinder on the batch of matches
    """
    result = []
    for match in matches:
        try:
            result.extend(run_single(**match, cached_compounds = cached_data, output_dir = output_dir))
        except Exception as err:
            result.append({"error": "Error in run_single: " + str(err), "match_index": match["match_index"]})
    
    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if save_pickle:
        with open(os.path.join(output_dir, file_name + ".pkl"), "wb") as f:
            pickle.dump(result, f)
    
    if save_csv:
        df = pd.DataFrame(result)
        df.to_csv(os.path.join(output_dir, file_name + ".csv"), index = False)
    
    return result


def runner(matches, output_dir = None, file_name = None, save_pickle = True, save_csv = True, cached_data = None, 
           number_of_cores = 1, batch_size = 1000, verbose = True) -> list:
    """
    run the modifinder on the batch of matches
    """
    
    if number_of_cores > 1:
        batches = list(chunkify(matches, batch_size))
        number_of_cores = min(number_of_cores, len(batches))
        
        # create a temporary directory to save the results
        if output_dir is None:
            output_dir = os.getcwd()
        temp_dir = os.path.join(output_dir, "temp")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers = number_of_cores) as executor:
            results = list(tqdm(executor.map(run_batch,
                                             batches,
                                             [temp_dir] * len(batches),
                                             [f"batch_{i}" for i in range(len(batches))],
                                             [save_pickle] * len(batches),
                                             [save_csv] * len(batches),
                                             [cached_data] * len(batches))
                                , total = len(batches), disable = not verbose))
            
        final_result = []
        for result in results:
            final_result.extend(result)
        
        # remove the temporary directory and its content
        for file in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, file))
        os.rmdir(temp_dir)
        
        
    else:
        final_result = run_batch(matches, output_dir, "batch", save_pickle = False, save_csv = False, cached_data = cached_data)
    
    
    if save_pickle:
        with open(os.path.join(output_dir, file_name + ".pkl"), "wb") as f:
            pickle.dump(final_result, f)
    
    if save_csv:
        df = pd.DataFrame(final_result)
        df.to_csv(os.path.join(output_dir, file_name + ".csv"), index = False)
        
    return final_result


def main(args: dict):
    """
    Run the modifinder algorithm on the given input data.
    """
    if args["input"].endswith(".csv"):
        data = pd.read_csv(args["input"])
        data = data.to_dict()
    elif args["input"].endswith(".pkl"):
        with open(args["input"], "rb") as f:
            data = pickle.load(f)
    elif args["input"].endswith(".json"):
        data = pd.read_json(args["input"])
        data = data.to_dict()
    else:
        raise ValueError("Please provide the input data in csv, pkl or json format.")
    
    if args["save_csv"] is False and args["save_pickle"] is False:
        raise ValueError("Please provide the format to save the output data.")
    
    if args["output"] is None:
        if args["verbose"]:
            print("Output directory is not provided. Saving the output in the current directory.")
        args["output"] = os.path.join(os.getcwd(), "output")
    
    if args["file_name"] is None:
        if args["verbose"]:
            print("Output file name is not provided. Saving the output as input file name.")
        args["file_name"] = os.path.basename(args["input"]).split(".")[0]
        
    
    if args["cached_data"] is not None:
        if args["cached_data"].endswith(".pkl"):
            with open(args["cached_data"], "rb") as f:
                cached_data = pickle.load(f)
        else:
            raise ValueError("Please provide the cached data in pkl format.")
    else:
        cached_data = None
    
    matches = []
    for i, match in enumerate(data):
        match["match_index"] = ("_" + str(i)) if "match_index" not in match else match["match_index"]
        matches.append(match)
    
    runner(matches, output_dir = args["output"], file_name = args["file_name"], save_pickle = args["save_pickle"],
              save_csv = args["save_csv"], cached_data = cached_data, number_of_cores = args["cores"],
              batch_size = args["batch_size"], verbose = args["verbose"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run the modifinder algorithm on the given input data.")
    parser.add_argument("input", help = "Path to the input data file.")
    parser.add_argument("--output", help = "Path to the output directory.")
    parser.add_argument("--file_name", help = "Name of the output file.")
    parser.add_argument("--cores", help = "Number of cores to use.", default = 1, type = int)
    parser.add_argument("--batch_size", help = "Size of the batch to run in parallel.", default = 1000, type = int)
    parser.add_argument("--verbose", help = "Print the progress bar.", dest="verbose", action="store_true")
    parser.add_argument("--no_verbose", help = "Do not print the progress bar.", dest="verbose", action="store_false")
    parser.add_argument("--save_pickle", help = "Save the output as a pickle file.", dest="save_pickle", action="store_true")
    parser.add_argument("--no_save_pickle", help = "Do not save the output as a pickle file.", dest="save_pickle", action="store_false")
    parser.add_argument("--save_csv", help = "Save the output as a csv file.", dest="save_csv", action="store_true")
    parser.add_argument("--no_save_csv", help = "Do not save the output as a csv file.", dest="save_csv", action="store_false")
    parser.add_argument("--cached_data", help = "Path to the cached data file.")
    
    parser.set_defaults(save_pickle = True, save_csv = True, verbose = False)
    args = parser.parse_args()
    
    args = vars(args)
    main(args)
